#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctime>

#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "tree.h"

#define BUFLEN 10000

using namespace std;

int main(int argc, char** argv) {
  #ifndef SKIMCHICO
  if(argc!=4) {
    cerr<<argv[0]<<" INPUTFILE OUTPUTFILE GTCalib"<<endl;
    exit(0);
  }
  #else
  if(argc!=5) {
    cerr<<argv[0]<<" INPUTFILE OUTPUTFILE GTCalib SKIMFILE"<<endl;
    exit(0);
  }
  #endif
//  string dir = "../../1563x_track/";
//  char filename[120];
//  sprintf(filename,"%sRun%04dTK.dat",dir.c_str(),atoi(argv[1]));
  FILE *in = fopen64(argv[1],"rb");
  if(!in){
    cerr<<" Cannot open file "<<argv[1]<<endl;
    return 1;
  }

  GEBheader aGeb;
  unsigned char cBuf[BUFLEN];
  
  int read;
  long long totread=0;
  int cnt=0;
  char option;
  long long int tsfirst=0;
  long long int tslast=0;

  TFile *opf=new TFile(argv[2],"RECREATE");
  TTree *opt=new TTree("t","opt");

  r2root *r2r=new r2root(opt);
  r2r->ReadGTCali(argv[3]);
  r2r->Clear();
  #ifdef SKIMCHICO
  r2r->OpenSkimFile(argv[4]);
  #endif

  opf->cd();
  clock_t starttime = clock();
  while(fread(&aGeb,sizeof(GEBheader),1,in)==1){
    totread+=sizeof(GEBheader);
    if(aGeb.length>BUFLEN){
      cerr<<endl;
      cerr<<"BUFLEN: "<<BUFLEN<<endl;
      cerr<<"GEBLEN: "<<aGeb.length<<endl;
      cerr<<"GEB type: "<<aGeb.type<<endl;
      exit(1);
    }
    if(aGeb.type != 3) aGeb.timestamp *= 10;
    if(tsfirst==0) {
      tsfirst = aGeb.timestamp;
    }
    if(tsfirst>aGeb.timestamp) tsfirst = aGeb.timestamp;
    if(tslast<aGeb.timestamp) tslast=aGeb.timestamp;
    if(TMath::Abs(tslast-tsfirst)>1000 * 10){
      if(r2r->filled()){
        r2r->CalTimeDiff();
        opt->Fill();
	#ifdef SKIMCHICO
	r2r->SkimChico();
	#endif
      }
      tsfirst=tslast=aGeb.timestamp;
      r2r->Clear();
    }
    read = fread(cBuf, sizeof(char), aGeb.length,in);
    totread+=read;
    if(aGeb.type != 3) aGeb.timestamp /= 10;
    if(aGeb.type==12){ // Chico
      r2r->ProcessChico(cBuf,&aGeb);
    }else if(aGeb.type== 1){ // Decomp
      r2r->ProcessGamma(cBuf,&aGeb);
    }else if(aGeb.type== 3){ // Track
      r2r->ProcessTKGamma(cBuf,&aGeb);
    }
    cnt++;
    if((cnt % 200)==0 ) {
      cerr << "Event " << cnt
	   << " read:" << read
	   << " total read:" << totread/1000000
	   << " Mb \r";
      cerr.flush();
    }
  }
  if(r2r->filled()){
    r2r->CalTimeDiff();
    opt->Fill(); // Fill possible last event
    #ifdef SKIMCHICO
    r2r->SkimChico();
    #endif
  }
  opt->Write();
  opf->Close();
  #ifdef SKIMCHICO
  r2r->CloseSkimFile();
  #endif
  starttime = clock()-starttime;
  cout<<endl;
  float durat = (float)starttime/CLOCKS_PER_SEC;
  cout<<"Used "<<durat<<" seconds"<<endl;
  cout<<"File size: "<<totread/1000000<<" Mb Avearge reading speed: "<< totread/1000000/(durat)<<" Mb/s"<<endl;
  cout<<"Converting "<<argv[1]<<" done!"<<endl;
}
