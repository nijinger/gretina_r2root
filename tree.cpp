#include <iostream>
#include <fstream>
#include <TVector3.h>
#include <TMath.h>
#include <TRandom.h>

#include "tree.h"

#define GTPosOffsetZ 0.0
#define GTPosOffsetY1 19.4
#define GTPosOffsetY2 19.4
#define PPAC_NUM 20

//#define CHICORAWDATA

using namespace std;

bool find_chico_full_singlehit(int anodech,int cnum,unsigned short ch[128], int& index){
  if(cnum<4||anodech*4<ch[0]||anodech*4>=ch[cnum-1]) return false;
  for(int i =0;i<cnum;i++){
    if(ch[i]%4!=0||anodech!=ch[i]/4) continue;
    else {
      if(i+4>cnum) return false;
      if(ch[i+1]-ch[i]!=1||ch[i+2]-ch[i+1]!=1||ch[i+3]-ch[i+2]!=1) return false;
      index=i;
      return true;
    }
  }
  return false;
}

bool find_chico_theta_singlehit(int anodech,int cnum,unsigned short ch[128], int& index){
  if(cnum<2||anodech*4<ch[0]||anodech*4>=ch[cnum-1]) return false;
  for(int i =0;i<cnum;i++){
    if(ch[i]%4!=0||anodech!=ch[i]/4) continue;
    else {
      if(i+2>cnum) return false;
      if(ch[i+1]-ch[i]!=1) return false;
      index=i;
      return true;
    }
  }
  return false;
}

bool recover_anode_by_cathode(int ach,int cnum, unsigned short int ch[128], int index, int &cindex, int& ocnum){
  if(cnum-index-4<2) return 0;
  int bch = (ach+5)%10;
  unsigned short och[128];
  for(int i=0;i<cnum-index-4;i++){
    och[i]=ch[i+index+4];
  }
  int tindex=-1;
  if(find_chico_full_singlehit(bch,cnum-index-4,och,tindex)){
    cindex = tindex+index+4;
    ocnum=4;
    return 1;
  }else if(find_chico_theta_singlehit(bch,cnum-index-4,och,tindex)){
    cindex = tindex+index+4;
    ocnum=2;
    return 1;
  }else return 0;
}

void r2root::chico_theta_phi(int ach, float *v, float &the, float &phi){
  the=gainTheta[ach]*(v[0]-v[1])+offTheta[ach]+gRandom->Uniform()-0.5;
  phi=gainPhi[ach]*(v[2]-v[3])+offPhi[ach]+gRandom->Uniform()-0.5;
  phi+=36.0*(float)(ach%10);
}


float crmat[MAXDETPOS+1][MAXCRYSTALNO+1][4][4];
void RdCRMAT(const char *fn) {

  FILE *fp;

  float f1, f2, f3, f4;
  int pos, xtal;
  int nn = 0;
  char *st, str[256];

  fp = fopen(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  }
  printf("\"%s\" open....", fn);

  /* Read values. */
  nn = 0;
  st = fgets(str, 256, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* Empty line, do nothing */
    } else {
      sscanf(str, "%i %i", &pos, &xtal);
      for (int i=0; i<4; i++) {
        st = fgets(str, 256, fp);
        sscanf(str, "%f %f %f %f", &f1, &f2, &f3, &f4);
        crmat[pos-1][xtal][i][0] = f1;
        crmat[pos-1][xtal][i][1] = f2;
        crmat[pos-1][xtal][i][2] = f3;
        crmat[pos-1][xtal][i][3] = f4;
      }
      nn++;
    }

    /* Attempt to read the next line. */
    st = fgets(str, 256, fp);
  }

  printf("Read %i rotation matrix coefficients.\n", nn);

  /* Done! */
  fclose(fp);
}

template<class T1,class T2>
void sort(int n,T1 val[128],T2 ch[128]){
  T1 tval[128];
  int index[128];
  T2 tch[128];
  memcpy(tval,val,sizeof(T1)*128);
  memcpy(tch,ch,sizeof(T2)*128);

  TMath::Sort(n,tch,index,0);

  for(int i=0;i<n;i++){
    val[i]=tval[index[i]];
    ch[i]=tch[index[i]];
  }
}

void ReadParticle(PARTICLE &particle,const BYTE* EventBuf,int length){
  const unsigned char *pos=EventBuf;
  int cnt=0;
  memcpy(&(particle.id),pos,4);
  pos+=4;
  cnt+=4;
  memcpy(&(particle.t),pos,56);
  pos+=56;
  cnt+=56;
  memcpy(&(particle.eR),pos,4);
  pos+=4;
  cnt+=4;
  memcpy(&(particle.rf),pos,8);
  pos+=8;
  cnt+=8;
  memcpy(&(particle.single),pos,1);
  pos+=4;
  cnt+=4;
  memcpy(&(particle.Mass),pos,4);
  pos+=4;
  cnt+=4;
  memcpy(&(particle.back),pos,1);
  cnt+=4;
  if((cnt!=length)){
    cout<<"cnt : "<<cnt<<endl;
    cout<<"length: "<<length<<endl;
  }
  assert(cnt==length);
}

r2root::r2root(TTree *iopt){
  assert(iopt!=0) ;
  opt=iopt;
  BranchOpt();
  RdCRMAT("crmat.dat");
//  InitPPACMap();
  ChicoInit();
}

void r2root::ChicoInit(){
  int i,j,k;
  float angle;
  string OneLine;

  ifstream THETACALFILE("ppacTheta.cal", ios::in);
  if(!THETACALFILE.is_open()) {
    cerr << "Error opening Calibration file ppacTheta.cal"<< endl;  
    exit(1);
  }
  getline(THETACALFILE,OneLine);
  //cout << OneLine << endl;
  getline(THETACALFILE,OneLine);
  //cout << OneLine << endl;
  for(i=0;i<PPAC_NUM;i++){
    THETACALFILE >> k >> offTheta[i] >> gainTheta[i];
    //THETACALFILE >> k >> offTheta[i] >> gainTheta[i] >> quadTheta[i];
    //cout << setw(3) << k << setw(15) << offTheta[i] << setw(15) << gainTheta[i] << endl;
  }
  THETACALFILE.close();

  ifstream PHICALFILE("ppacPhi.cal", ios::in);
  if(!PHICALFILE.is_open()) {
    cerr << "Error opening Calibration file ppacPhi.cal "<< endl;  
    exit(1);
  }
  getline(PHICALFILE,OneLine);
  //cout << OneLine << endl;
  getline(PHICALFILE,OneLine);
  //cout << OneLine << endl;
  for(i=0;i<PPAC_NUM;i++){
    PHICALFILE >> k >> offPhi[i] >> gainPhi[i];
    //cout << setw(3) << k << setw(15) << offPhi[i] << setw(15) << gainPhi[i] << endl;
  }
  PHICALFILE.close();

}

bool r2root::BranchOpt(){

//  opt->Branch("gev.",&gev);

  opt->Branch("np",&np,"np/I");
  opt->Branch("cts",&cts,"cts[np]/l");
  opt->Branch("hts",&hts,"hts[np]/l");
  opt->Branch("id",&id,"id[np]/I");
  opt->Branch("dT",&dT,"dT[np]/I");
  opt->Branch("dL",&dL,"dL[np]/F");
  opt->Branch("dR",&dR,"dR[np]/F");
  opt->Branch("phiL",&phiL,"phiL[np]/I");
  opt->Branch("phiR",&phiR,"phiR[np]/I");
  opt->Branch("fphiL",&fphiL,"fphiL[np]/F");
  opt->Branch("fphiR",&fphiR,"fphiR[np]/F");
  opt->Branch("thetaL",&thetaL,"thetaL[np]/I");
  opt->Branch("thetaR",&thetaR,"thetaR[np]/I");
  opt->Branch("fthetaL",&fthetaL,"fthetaL[np]/F");
  opt->Branch("fthetaR",&fthetaR,"fthetaR[np]/F");

  opt->Branch("ng",&ng,"ng/I");
//  opt->Branch("gts",&gts,"gts[ng]/L");
  opt->Branch("x",&x,"x[ng]/F");
  opt->Branch("y",&y,"y[ng]/F");
  opt->Branch("z",&z,"z[ng]/F");
  opt->Branch("e",&e,"e[ng]/F");
  opt->Branch("theta",&theta,"theta[ng]/F");
  opt->Branch("phi",&phi,"phi[ng]/F");
  opt->Branch("t",&t,"t[ng]/l");
  opt->Branch("offt",&offt,"offt[ng]/D");
  opt->Branch("t0",&t0,"t0[ng]/D");
  opt->Branch("dtpg",&dtpg,"dtpg[ng]/D");
  opt->Branch("cc_id",&cc_id,"cc_id[ng]/I");
  opt->Branch("cry_id",&cry_id,"cry_id[ng]/I");

  opt->Branch("ntg",&ntg,"ntg/I");
  opt->Branch("tsTK",&tsTK,"tsTK/l");
  opt->Branch("pad",&pad,"pad[ntg]/I");
  opt->Branch("tracked",&tracked,"tracked[ntg]/I");
  opt->Branch("esum",&esum,"esum[ntg]/F");
  opt->Branch("ndet",&ndet,"ndet[ntg]/I");
  opt->Branch("gtkts",&gtkts,"gtkts[ntg]/l");
  opt->Branch("fom",&fom,"fom[ntg]/F");
  opt->Branch("x0",&x0,"x0[ntg]/F");
  opt->Branch("y0",&y0,"y0[ntg]/F");
  opt->Branch("z0",&z0,"z0[ntg]/F");
  opt->Branch("e0",&e0,"e0[ntg]/F");
  opt->Branch("x1",&x1,"x1[ntg]/F");
  opt->Branch("y1",&y1,"y1[ntg]/F");
  opt->Branch("z1",&z1,"z1[ntg]/F");
  opt->Branch("e1",&e1,"e1[ntg]/F");
  opt->Branch("fhcrID",&fhcrID,"fhcrID[ntg]/I");
  opt->Branch("gid",&gid,"gid[ntg]/I");
  opt->Branch("gofft",&gofft,"gofft[ntg]/D");
  opt->Branch("dtptg",&dtptg,"dtptg[ntg]/D");
/*
  opt->Branch("chico.LEDts",&ChicoEvent,"LEDts/l");
  opt->Branch("chico.cathode_tdc_val",&ChicoEvent.cathode_tdc_val[0],"cathode_tdc_val[128]/I");
  opt->Branch("chico.anode_tdc_val",&ChicoEvent.anode_tdc_val[0],"anode_tdc_val[128]/I");
  opt->Branch("chico.anode_qdc_val",&ChicoEvent.anode_qdc_val[0],"anode_qdc_val[128]/I");
  opt->Branch("chico.status",&ChicoEvent.status,"status/i");
  opt->Branch("chico.RF",&ChicoEvent.RF,"RF/I");
  opt->Branch("chico.cathode_tdc_num",&ChicoEvent.cathode_tdc_num,"cathode_tdc_num/s");
  opt->Branch("chico.cathode_tdc_ch",&ChicoEvent.cathode_tdc_ch[0],"cathode_tdc_ch[128]/s");
  opt->Branch("chico.anode_tdc_num",&ChicoEvent.anode_tdc_num,"anode_tdc_num/s");
  opt->Branch("chico.anode_tdc_ch",&ChicoEvent.anode_tdc_ch[0],"anode_tdc_ch[128]/s");
  opt->Branch("chico.anode_qdc_num",&ChicoEvent.anode_qdc_num,"anode_qdc_num/s");
  opt->Branch("chico.anode_qdc_ch",&ChicoEvent.anode_qdc_ch[0],"anode_qdc_ch[128]/s");
  opt->Branch("chico.SINGLE",&ChicoEvent.SINGLE,"SINGLE/O");
*/

#ifdef CHICORAWDATA
  opt->Branch("chico.LEDts",&ChicoEvent,"LEDts/l");
  opt->Branch("chico.ctv",&ChicoEvent.cathode_tdc_val[0],"cathode_tdc_val[128]/I");
  opt->Branch("chico.atv",&ChicoEvent.anode_tdc_val[0],"anode_tdc_val[128]/I");
  opt->Branch("chico.aqv",&ChicoEvent.anode_qdc_val[0],"anode_qdc_val[128]/I");
  opt->Branch("chico.status",&ChicoEvent.status,"status/i");
  opt->Branch("chico.RF",&ChicoEvent.RF,"RF/I");
  opt->Branch("chico.ctn",&ChicoEvent.cathode_tdc_num,"cathode_tdc_num/s");
  opt->Branch("chico.ctc",&ChicoEvent.cathode_tdc_ch[0],"cathode_tdc_ch[128]/s");
  opt->Branch("chico.atn",&ChicoEvent.anode_tdc_num,"anode_tdc_num/s");
  opt->Branch("chico.atc",&ChicoEvent.anode_tdc_ch[0],"anode_tdc_ch[128]/s");
  opt->Branch("chico.aqn",&ChicoEvent.anode_qdc_num,"anode_qdc_num/s");
  opt->Branch("chico.aqc",&ChicoEvent.anode_qdc_ch[0],"anode_qdc_ch[128]/s");
  opt->Branch("chico.SINGLE",&ChicoEvent.SINGLE,"SINGLE/O");
#endif

  opt->Branch("nhit",&nhit,"nhit/I");
  opt->Branch("ach",ach,"ach[nhit]/I");
  opt->Branch("aval",aval,"aval[nhit]/F");
  opt->Branch("chit",&chit,"chit/I");
  opt->Branch("cch",cch,"cch[chit]/I");
  opt->Branch("cval",cval,"cval[chit]/F");
/*
  opt->Branch("tthe",tthe,"tthe[2]");
  opt->Branch("tphi",tphi,"tphi[2]");
*/
}

void r2root::Clear(){
  np=0;
  ng=0;
  ntg=0;
  tsTK=-1;
//  memset(&gev,0,sizeof(gev));
  memset(&ChicoEvent,0,sizeof(ChicoEvent));
  memset(tthe,0,sizeof(tthe));
  memset(tphi,0,sizeof(tphi));
  nhit=0;
  chit=0;
  #ifdef SKIMCHICO
  memset(&Skim,0,sizeof(Skim));
  #endif
}

int r2root::GetChicoEvent(const unsigned int *DataRecord, CHICOEVENT *ChicoEvent){
  unsigned short int EvSize;
  unsigned short int chan=0;
  int val=0,refval=0;
  unsigned int NextInt;
  int seenTrailer=0;
  int i,j,k=0;
  int CH_Counter=0;
  static int multiAnodeTDCNum=0;
  static int multiCathodeTDCNum=0;

  EvSize = DataRecord[0]/4;
  ChicoEvent->status = (unsigned short int)((DataRecord[0] & 0xffff0000)>>16);
  ChicoEvent->LEDts = (unsigned long long int) (DataRecord[1] & 0xffffffff);
  ChicoEvent->LEDts += ((unsigned long long int) (DataRecord[2] & 0xffff)<<32);
  j = 2;
  NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
  j++;
  if(NextInt != 0xffffffff){					//Anode_QDC
    //assert(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) == QDCHEADER);
    if(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT)!=QDCHEADER)return 0;
    //assert(((NextInt & QDCGEOMASK) >> QDCGEOSHIFT) == ANODE_E_VSN);
    if(((NextInt & QDCGEOMASK) >> QDCGEOSHIFT) != ANODE_E_VSN)return 0;
    CH_Counter = (NextInt & COUNTMASK) >> COUNTSHIFT;
    k=0;
    for(i=0;i<CH_Counter;i++){
      NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
      j++;
      //assert(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) == DATA);
      if(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) != DATA)return 0;
      chan = (unsigned short int)((NextInt & QDCCHANMASK) >> QDCCHANSHIFT);
      val = (NextInt & QDCDATAMASK);
      if(chan<PPAC_NUM && val >0 && k<32){
        ChicoEvent->anode_qdc_ch[k] = chan;
        ChicoEvent->anode_qdc_val[k] = val;
        k++;
      }
    }
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    //assert(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) == QDCTRAILER);
    if(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) != QDCTRAILER)return 0;
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    if(NextInt != 0xffffffff){
      while (1){
        NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
        j++;
        if(NextInt == 0xffffffff)break;
      }
    }
  }
  ChicoEvent->anode_qdc_num = k;
  CH_Counter=0;
  chan=0;val=0;
  ChicoEvent->SINGLE = false;
  NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
  j++;
  if(NextInt != 0xffffffff){					//Anode_TDC
    //assert((NextInt & TDCTYPEMASK) == TDCHEADER);
    if((NextInt & TDCTYPEMASK) != TDCHEADER)return 0;
    //assert((NextInt & TDCGEOMASK) == ANODE_T_VSN);
    if((NextInt & TDCGEOMASK) != ANODE_T_VSN)return 0;
    while(!seenTrailer){
      NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
      j++;
      switch(NextInt & TDCTYPEMASK){
        case DATA:
          chan = (unsigned short int)((NextInt & TDCCHANMASK) >> TDCCHANSHIFT); 
          val =(NextInt & TDCDATAMASK);
          if(chan != ANODE_REFCH && chan != RFCH && chan != SingleFlag){
            if(chan<PPAC_NUM && CH_Counter<128){
              ChicoEvent->anode_tdc_ch[CH_Counter] = chan;
              ChicoEvent->anode_tdc_val[CH_Counter] =val;
              CH_Counter++;
            }
          }
          else if (chan == RFCH){
            ChicoEvent->RF = val;
          }
          else if (chan == SingleFlag){
            ChicoEvent->SINGLE = true;
          }
          else if(chan == ANODE_REFCH){
            refval = (NextInt & TDCDATAMASK);
          }
          break;
        case TDCTRAILER:
          seenTrailer = 1;
          break;
        default:
          break;
      }
    }
    //if(refval > 0){
      for(i=0;i<CH_Counter;i++){
        ChicoEvent->anode_tdc_val[i] -= refval;
      }
      //ChicoEvent->RF -= refval;
    //}
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    if(NextInt != 0xffffffff) {
      multiAnodeTDCNum++;
      if((multiAnodeTDCNum%100000)==0)printf("*Warning* anode TDC package with multiple events: %i\n",multiAnodeTDCNum);
      while (1){
        NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
        j++;
        if(NextInt == 0xffffffff)break;
      }
    }
  }
  ChicoEvent->anode_tdc_num = CH_Counter;
  CH_Counter=0; seenTrailer=0;
  chan=0;val=0;
  NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
  j++;
  if(NextInt != 0xffffffff){					//Cathode_TDC
    //assert((NextInt & TDCTYPEMASK) == TDCHEADER);
    if((NextInt & TDCTYPEMASK) != TDCHEADER)return 0;
    //assert((NextInt & TDCGEOMASK) == CATHODE_T_VSN);
    if((NextInt & TDCGEOMASK) != CATHODE_T_VSN)return 0;
    while(!seenTrailer){
      NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
      j++;
      switch(NextInt & TDCTYPEMASK){
        case DATA:
          chan = (unsigned short int)((NextInt & TDCCHANMASK) >> TDCCHANSHIFT); 
          val =(NextInt & TDCDATAMASK);
          if(chan != CATHODE_REFCH){
            if(chan < PPAC_NUM*4 && CH_Counter<128){
	      if(chan==57) chan=59;
	      else if(chan==59) chan=57;
              ChicoEvent->cathode_tdc_ch[CH_Counter] = chan;
              ChicoEvent->cathode_tdc_val[CH_Counter] = val;
              CH_Counter++;
            }
          }
          else if(chan == CATHODE_REFCH){
            refval = (NextInt & TDCDATAMASK);
          }
          break;
        case TDCTRAILER:
          seenTrailer = 1;
          break;
        default:
          break;
      }
    }
    //if(refval > 0){
      //for(i=0;i<CH_Counter;i++){
      //  ChicoEvent->cathode_tdc_val[i] -= refval;
      //}
    //}
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    if(NextInt != 0xffffffff) {
      multiCathodeTDCNum++;
      if((multiCathodeTDCNum%100000)==0)printf("*Warning* cathode TDC package with multiple events: %i\n",multiCathodeTDCNum);
      while (1){
        NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
        j++;
        if(NextInt == 0xffffffff)break;
      }
    }
  }
  ChicoEvent->cathode_tdc_num = CH_Counter;
  return 1;
}

void r2root::ProcessChico(const BYTE *buf,GEBheader *header){
   if(header->length>=1024) return;
   GetChicoEvent((const unsigned int*)buf,&ChicoEvent);
   sort(ChicoEvent.anode_tdc_num,ChicoEvent.anode_tdc_val,ChicoEvent.anode_tdc_ch);
   sort(ChicoEvent.cathode_tdc_num,ChicoEvent.cathode_tdc_val,ChicoEvent.cathode_tdc_ch);

   PARTICLE chico;
//   GetParticle(&ChicoEvent,&chico);
   if(GetParticle(ChicoEvent,chico)==0) return;
//   ReadParticle(chico,buf,header->length);
   if(np<=0) np=0;
   else if (np>=MaxChicoNum) return;
   hts[np]=header->timestamp;
   cts[np]=chico.t;
   id[np]=chico.id;
   dT[np]=chico.dT;
   dL[np]=chico.dL;
   dR[np]=chico.dR;
   phiL[np]=chico.phiL;
   phiR[np]=chico.phiR;
   fphiL[np]=chico.fphiL;
   fphiR[np]=chico.fphiR;
   thetaL[np]=chico.thetaL;
   thetaR[np]=chico.thetaR;
   fthetaL[np]=chico.fthetaL;
   fthetaR[np]=chico.fthetaR;
   #ifdef SKIMCHICO
//   memcpy(&Chico[np].header,header,sizeof(GEBheader));
   Chico[np].header.type=header->type;
   Chico[np].header.length=sizeof(PARTICLE);
   Chico[np].header.timestamp=header->timestamp;
   Chico[np].particle=chico;
   #endif
   np++;
}


int r2root::GetParticle(CHICOEVENT *ChicoEvent, PARTICLE *particle){
}

int r2root::GetParticle(CHICOEVENT &ChicoEvent, PARTICLE &particle){
  
  //No valid anode information was recorded
  if(ChicoEvent.anode_tdc_num==0||ChicoEvent.cathode_tdc_num<4) return 0;

  int cnum=0;
  for(int i=0;i<ChicoEvent.anode_tdc_num;i++){
    int index = 0;
    if(find_chico_full_singlehit(ChicoEvent.anode_tdc_ch[i],ChicoEvent.cathode_tdc_num,ChicoEvent.cathode_tdc_ch,index)){
      ach[nhit]=ChicoEvent.anode_tdc_ch[i];
      aval[nhit]=ChicoEvent.anode_tdc_val[i];
      cch[nhit*4]=ChicoEvent.cathode_tdc_ch[index];
      cch[nhit*4+1]=ChicoEvent.cathode_tdc_ch[index+1];
      cch[nhit*4+2]=ChicoEvent.cathode_tdc_ch[index+2];
      cch[nhit*4+3]=ChicoEvent.cathode_tdc_ch[index+3];
      cval[nhit*4]=ChicoEvent.cathode_tdc_val[index];
      cval[nhit*4+1]=ChicoEvent.cathode_tdc_val[index+1];
      cval[nhit*4+2]=ChicoEvent.cathode_tdc_val[index+2];
      cval[nhit*4+3]=ChicoEvent.cathode_tdc_val[index+3];
      nhit++;
      cnum+=4;
    }else if(find_chico_theta_singlehit(ChicoEvent.anode_tdc_ch[i],ChicoEvent.cathode_tdc_num,ChicoEvent.cathode_tdc_ch,index)){
      ach[nhit]=ChicoEvent.anode_tdc_ch[i];
      aval[nhit]=ChicoEvent.anode_tdc_val[i];
      cch[nhit*4]=ChicoEvent.cathode_tdc_ch[index];
      cch[nhit*4+1]=ChicoEvent.cathode_tdc_ch[index+1];
      cch[nhit*4+2]=ChicoEvent.cathode_tdc_ch[index]+2;
      cch[nhit*4+3]=ChicoEvent.cathode_tdc_ch[index]+3;
      cval[nhit*4]=ChicoEvent.cathode_tdc_val[index];
      cval[nhit*4+1]=ChicoEvent.cathode_tdc_val[index+1];
      cval[nhit*4+2]=-999;
      cval[nhit*4+3]=-999;
      nhit++;
      cnum+=2;
    }
  }

  // if more than 2 remaning catodes events have not been correctly assigned, try to figure out missing anode events
  if(ChicoEvent.cathode_tdc_num-cnum>2){
    for(int i = 0;i<ChicoEvent.cathode_tdc_num-1;i++){
      int index = -1;
      if(ChicoEvent.cathode_tdc_ch[i]%4==0){
        int flag=0;
	int tach = ChicoEvent.cathode_tdc_ch[i]/4;
        for(int j=0;j<nhit;j++){
	  if(tach==ach[j]) flag=1;
	}
	if(flag==0&&ChicoEvent.cathode_tdc_ch[i+1]-ChicoEvent.cathode_tdc_ch[i]==1){
	  if(find_chico_full_singlehit(tach,ChicoEvent.cathode_tdc_num,ChicoEvent.cathode_tdc_ch,index)){
            ach[nhit]=tach;
            aval[nhit]=-999;
            cch[nhit*4]=ChicoEvent.cathode_tdc_ch[index];
            cch[nhit*4+1]=ChicoEvent.cathode_tdc_ch[index+1];
            cch[nhit*4+2]=ChicoEvent.cathode_tdc_ch[index+2];
            cch[nhit*4+3]=ChicoEvent.cathode_tdc_ch[index+3];
            cval[nhit*4]=ChicoEvent.cathode_tdc_val[index];
            cval[nhit*4+1]=ChicoEvent.cathode_tdc_val[index+1];
            cval[nhit*4+2]=ChicoEvent.cathode_tdc_val[index+2];
            cval[nhit*4+3]=ChicoEvent.cathode_tdc_val[index+3];
            nhit++;
            cnum+=4;
	  }else if(find_chico_theta_singlehit(tach,ChicoEvent.cathode_tdc_num,ChicoEvent.cathode_tdc_ch,index)){
            ach[nhit]=tach;
            aval[nhit]=-999;
            cch[nhit*4]=ChicoEvent.cathode_tdc_ch[index];
            cch[nhit*4+1]=ChicoEvent.cathode_tdc_ch[index+1];
            cch[nhit*4+2]=ChicoEvent.cathode_tdc_ch[index]+2;
            cch[nhit*4+3]=ChicoEvent.cathode_tdc_ch[index]+3;
            cval[nhit*4]=ChicoEvent.cathode_tdc_val[index];
            cval[nhit*4+1]=ChicoEvent.cathode_tdc_val[index+1];
            cval[nhit*4+2]=-999;
            cval[nhit*4+3]=-999;
            nhit++;
            cnum+=2;
	  }
	}
      }

    }
  }

  if(nhit==0) return 0;
  chit=nhit*4;

  sort(nhit,aval,ach);
  sort(nhit*4,cval,cch);

  int flag = 0;
  for(int i=0;i<nhit;i++){
    if(cval[i*4+2]!=-999) {
      flag=1;
      break;
    }
  }
  //cannot find one fully recovered event
  if(flag==0) return 0;

  // only found one single front hit
  if(nhit==1&&ach[0]<10) return 0;
  
  float theta[2],phi[2];
  int index1=-1,index2=-1;
  // one single back hit
  if(nhit==1&&ach[0]>9){
    chico_theta_phi(ach[0],&cval[0],theta[0],phi[0]);
    tthe[0]=theta[0];
    tphi[0]=phi[0];
  //  return 1;
  }else if(nhit>=2){
    for(int i=0;i<nhit-1;i++)
    for(int j=i+1;j<nhit;j++){
      if(ach[j]-ach[i]==5||ach[j]-ach[i]==15) {
        index1=i;index2=j;
      }
    }
    if(index1==-1||(cval[index1*4+2]==-999&&cval[index2*4+2]==-999)) return 0;
    if(ach[index2]%10<ach[index1]%10) swap(index1,index2);
    chico_theta_phi(ach[index1],&cval[index1*4],theta[0],phi[0]);
    chico_theta_phi(ach[index2],&cval[index2*4],theta[1],phi[1]);
    if(cval[index1*4+2]==-999){
      phi[0]=phi[1]-180.;
    }else if(cval[index2*4+2]==-999){
      phi[1]=phi[0]+180.;
//	if(phi[1]>360) phi[1]-=360.;
    }
    tthe[0]=theta[0];
    tphi[0]=phi[0];
    tthe[1]=theta[1];
    tphi[1]=phi[1];
//    return 1;
  }
  
  if(nhit==1){
     particle.t= (double)ChicoEvent.LEDts;
     particle.dT=-999;
     particle.id=ach[0];
     particle.fthetaL=theta[0]*M_PI/180.;
     particle.fphiL=phi[0]*M_PI/180.;
     particle.single = ChicoEvent.SINGLE;
     particle.back = true;
     return 1;
  }

  if(aval[index1]!=-999&&aval[index2]!=-999){
    particle.t= (double)ChicoEvent.LEDts;
    particle.rf = ((double)ChicoEvent.RF * ch2ns)*0.1 + gRandom->Uniform()-0.5;
    particle.id=ach[index1];
    particle.dT=(aval[index1]-aval[index2]);
    particle.fthetaL=theta[0]*M_PI/180.;
    particle.fthetaR=theta[1]*M_PI/180.;
    particle.fphiL = phi[0]*M_PI/180.;
    particle.fphiR = phi[1]*M_PI/180.;
    float dL = 12.8/(0.75471*sinf(particle.fthetaL)*cosf(particle.fphiL-(18.0+(float)(ach[index1]%10)*36.0)*M_PI/180.)
        +0.65606*cosf(particle.fthetaL));
    float dR = 12.8/(0.75471*sinf(particle.fthetaR)*cosf(particle.fphiR-(18.0+(float)(ach[index2]%10)*36.0)*M_PI/180.)
        +0.65606*cosf(particle.fthetaR));
    particle.dL = dL;
    particle.dR = dR;
    particle.single = ChicoEvent.SINGLE;
    particle.back=false;
    return 1;
  }else if(ach[index2]>9||ach[index1]>9){
     particle.t= (double)ChicoEvent.LEDts;
     particle.dT=-999;
     particle.id=ChicoEvent.anode_tdc_ch[index1];
     particle.fthetaL=theta[0]*M_PI/180.;
     particle.fthetaR=theta[1]*M_PI/180.;
     particle.fphiL = phi[0]*M_PI/180.;
     particle.fphiR = phi[1]*M_PI/180.;
     particle.single = ChicoEvent.SINGLE;
     particle.back=true;
     return 1;
  }
  return 0;
}
void r2root::ProcessTKGamma(const BYTE *buf,GEBheader *header){
   TRACKED_GAMMA_HIT gev;
   tsTK = header->timestamp;
   memcpy(&gev,buf,header->length);
   if(ntg<0) ntg=0;
   else if(ntg>=MaxGTNum) return;
   for(int i=0;i<gev.ngam;i++){
//     if(gev.pad!=0) continue;
     pad[ntg]=gev.pad;
     esum[ntg]=gev.gr[i].esum;
     ndet[ntg]=gev.gr[i].ndet;
     fom[ntg]=gev.gr[i].fom;
     tracked[ntg]=gev.gr[i].tracked;
     gtkts[ntg]=gev.gr[i].timestamp;
     x0[ntg]=gev.gr[i].x0;
     y0[ntg]=gev.gr[i].y0;
     z0[ntg]=gev.gr[i].z0;
     e0[ntg]=gev.gr[i].e0;
     x1[ntg]=gev.gr[i].x1;
     y1[ntg]=gev.gr[i].y1;
     z1[ntg]=gev.gr[i].z1;
     e1[ntg]=gev.gr[i].e1;
#if (TRACK2 == 1)
     fhcrID[ntg] = gev.gr[i].fhcrID;
     int detectorPosition = ((fhcrID[ntg] & 0xfffc)>>2)-1;
     int crystalNumber = (fhcrID[ntg] & 0x0003);
     int cc_id = crystalNumber;
     int cry_id = -1;
     switch(detectorPosition+1){
       case Q1:
         cry_id=1;
         break;
       case Q2:
         cry_id=2;
         break;
       case Q3:
         cry_id=3;
         break;
       case Q4:
         cry_id=4;
         break;
       case Q5:
         cry_id=5;
         break;
       case Q6:
         cry_id=6;
         break;
       case Q7:
         cry_id=7;
         break;
       case Q8:
         cry_id=8;
         break;
       case Q9:
         cry_id=9;
         break;
       case Q10:
         cry_id=10;
         break;
       case Q11:
         cry_id=11;
         break;
       default:
         break;
     }

     gid[ntg] = (cry_id-1)*4+cc_id;
     gofft[ntg] = offdt[gid[ntg]];
#endif
     ntg++;
     if(ntg>=MaxGTNum) return;
   }
}

void r2root::ProcessGamma(const BYTE *buf,GEBheader *header){
  GTEVENT2 gev;
  assert(header->length == sizeof(GTEVENT2));
  memcpy(&gev,buf,header->length);
//  if(gev.pad!=0) return; // unsucessful decomposition
  if(ng<0) ng=0;
  else if(ng>=MaxGTNum) return;
  #ifdef SKIMCHICO
  memcpy(&Mode2[ng].header,header,sizeof(GEBheader));
//  cout<<"Mode2: "<<Mode2[ng].header.length<<endl;
  Mode2[ng].event2=gev;
  #endif
//  gts[ng]=header->timestamp;
  int nmax=0;
  float emax=0;
  for(int i=0;i<gev.num;i++)
  if(gev.intpts[i].e>emax){
    emax=gev.intpts[i].e;
    nmax=i;
  }
  int detectorPosition = ((gev.crystal_id & 0xfffc)>>2)-1;
  int crystalNumber = (gev.crystal_id & 0x0003);
  cc_id[ng] = crystalNumber;
  switch(detectorPosition+1){
    case Q1:
      cry_id[ng]=1;
      break;
    case Q2:
      cry_id[ng]=2;
      break;
    case Q3:
      cry_id[ng]=3;
      break;
    case Q4:
      cry_id[ng]=4;
      break;
    case Q5:
      cry_id[ng]=5;
      break;
    case Q6:
      cry_id[ng]=6;
      break;
    case Q7:
      cry_id[ng]=7;
      break;
    case Q8:
      cry_id[ng]=8;
      break;
    case Q9:
      cry_id[ng]=9;
      break;
    case Q10:
      cry_id[ng]=10;
      break;
    case Q11:
      cry_id[ng]=11;
 default:
      break;
  }
  x[ng] = ( (crmat[detectorPosition][crystalNumber][0][0] * gev.intpts[nmax].x) +
               (crmat[detectorPosition][crystalNumber][0][1] * gev.intpts[nmax].y) +
               (crmat[detectorPosition][crystalNumber][0][2] * gev.intpts[nmax].z) +
               (crmat[detectorPosition][crystalNumber][0][3]) );

  y[ng] = ( (crmat[detectorPosition][crystalNumber][1][0] * gev.intpts[nmax].x) +
               (crmat[detectorPosition][crystalNumber][1][1] * gev.intpts[nmax].y) +
               (crmat[detectorPosition][crystalNumber][1][2] * gev.intpts[nmax].z) +
               (crmat[detectorPosition][crystalNumber][1][3]) );

  z[ng] = ( (crmat[detectorPosition][crystalNumber][2][0] * gev.intpts[nmax].x) +
               (crmat[detectorPosition][crystalNumber][2][1] * gev.intpts[nmax].y) +
               (crmat[detectorPosition][crystalNumber][2][2] * gev.intpts[nmax].z) +
               (crmat[detectorPosition][crystalNumber][2][3]) );

  TVector3 dir(x[ng],y[ng],z[ng]);
  theta[ng]=dir.Theta();
  phi[ng]=dir.Phi();
  
  int id = (cry_id[ng]-1)*4+cc_id[ng];
  e[ng] = gev.tot_e*gainGT[id]+offGT[id];

  //gamma->t = (double)GTEvent2->timestamp + (double)GTEvent2->t0/10.;
  t[ng] = gev.timestamp;
  offt[ng] = (double)offdt[id];
  t0[ng] = (double)gev.t0;
  
  ng++;
}

void r2root::ReadGTCali(const char *GTCalib){
  ifstream GTCALFILE(GTCalib, ios::in);
  if(!GTCALFILE.is_open()) {
    cerr << "Error opening Gretina Calibration file!"<< endl;  
    exit(1);
  }
  for(int i=0;i<QNum*4;i++){
    GTCALFILE >> offGT[i] >> gainGT[i] >> offGT2[i] >> gainGT2[i] >> offdt[i];
    //cout << setw(3) << i+1 << setw(15) << offGT[i] << setw(15) << gainGT[i] << setw(15) << offdt[i] << endl;
  }
  GTCALFILE.close();
}

void r2root::CalTimeDiff(){
  for(int i = 0;i<ng;i++){
    double tmpt = (double)cts[0]-(double)t[i];
    dtpg[i] = tmpt-t0[i]-offt[i];
  }
  for(int i = 0;i<ntg;i++){
    double tmpt = (double)cts[0]*10-(double)gtkts[i];
    dtptg[i] = tmpt-gofft[i]*10;
  }
}

#ifdef SKIMCHICO

void r2root::OpenSkimFile(const char* skimf){
  fskim = fopen(skimf,"wb");
  if(fskim==NULL){
    cout<<"Cannot open skim file!"<<endl;
    exit(0);
  }
//  cout<<"fskim:"<<fskim<<endl;
  TFile *topf = new TFile("etCut.root");
  etCut = (TCutG*)topf->Get("etCut");
  topf->Close();
}

void r2root::CloseSkimFile(){
  fclose(fskim);
}

void r2root::SkimChico(){
  if(np<=0||ng<=0) return;
  for(int i=0;i<np;i++){
    Skim.chico[Skim.nCoinChico]=Chico[i];
    Skim.nCoinChico++;
  }
  for(int i=0;i<ng;i++){
//    if(abs(dtpg[i])<5){
    if(etCut->IsInside(dtpg[i],e[i])){
      Skim.mode2[Skim.nCoinGT]=Mode2[i];
      Skim.nCoinGT++;
    }
  }
  if(Skim.nCoinChico<=0||Skim.nCoinGT<=0) return;

  for(int i=0;i<Skim.nCoinGT;i++){
    fwrite(&Skim.mode2[i],sizeof(MODE2),1,fskim);
  }
  for(int i=0;i<Skim.nCoinChico;i++){
    fwrite(&Skim.chico[i],sizeof(CHICO),1,fskim);
  }
  
}
#endif
