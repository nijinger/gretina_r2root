#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <assert.h>
#include "GEBSort.h"

#define PPAC_NUM 20
#define ch2ns 0.100

//#define SKIMCHICO

typedef unsigned char BYTE;

class r2root {
public:
  r2root(){};
  r2root(TTree *opt);
  virtual ~r2root(){}

  void ProcessChico(const BYTE *buf, GEBheader *header);
  void ProcessTKGamma(const BYTE *buf,GEBheader *header);
  void ProcessGamma(const BYTE *buf,GEBheader *header);

  void Clear();
  bool BranchOpt();
  bool filled(){return (np>0||ng>0||ntg>0);}
  int GetChicoEvent(const unsigned int*,CHICOEVENT*);
  int GetParticle(CHICOEVENT *ChicoEvent, PARTICLE* particle);
  int GetParticle(CHICOEVENT &ChicoEvent, PARTICLE& particle);
  void ReadGTCali(const char *);
  void CalTimeDiff();
#ifdef SKIMCHICO
  void OpenSkimFile(const char*);
  void SkimChico();
  void CloseSkimFile();
#endif


private:
  void ChicoInit();
  void chico_theta_phi(int ach, float *v, float &the, float &phi);

  TTree *opt;

  Int_t np; // particle number
  ULong64_t cts[MaxChicoNum]; // chico ts from GEB
  ULong64_t hts[MaxChicoNum];
  Int_t id[MaxChicoNum]; // id of ppac
  Int_t dT[MaxChicoNum]; // time of flight difference
  Float_t dL[MaxChicoNum]; // distance to left
  Float_t dR[MaxChicoNum]; // distance to right
  Int_t phiL[MaxChicoNum]; // phi angle of left
  Int_t phiR[MaxChicoNum]; // phi angle of right
  Float_t fphiL[MaxChicoNum];
  Float_t fphiR[MaxChicoNum];
  Int_t thetaL[MaxChicoNum];
  Int_t thetaR[MaxChicoNum];
  Float_t fthetaL[MaxChicoNum];
  Float_t fthetaR[MaxChicoNum];

  Int_t ng; // gamma number
//  Long64_t gts[10]; // timestamp from GEB, identical with below t0, removed from entry
  Float_t x[MaxGTNum],y[MaxGTNum],z[MaxGTNum],e[MaxGTNum];
  Float_t theta[MaxGTNum],phi[MaxGTNum];
  ULong64_t t[MaxGTNum];
  Double_t t0[MaxGTNum],offt[MaxGTNum],dtpg[MaxGTNum];
  Int_t cc_id[MaxGTNum],cry_id[MaxGTNum];

  Int_t ntg;// traked gamma number
  ULong64_t tsTK; // global timestamp for tracked gamma
  Int_t pad[MaxGTNum]; // decomposition indicator
  Int_t tracked[MaxGTNum]; // 1 means sucessfully tracked
  Float_t esum[MaxGTNum]; // gamma ray energy
  Int_t ndet[MaxGTNum]; // number of interaction points
  ULong64_t gtkts[MaxGTNum]; // timestap of 1st interaction point
  Float_t fom[MaxGTNum]; // figure of merit
  Float_t x0[MaxGTNum],y0[MaxGTNum],z0[MaxGTNum],e0[MaxGTNum]; // 1st inter. point
  Float_t x1[MaxGTNum],y1[MaxGTNum],z1[MaxGTNum],e1[MaxGTNum]; // 2nd inter. point
  Int_t fhcrID[MaxGTNum];
  Int_t gid[MaxGTNum];
  Double_t gofft[MaxGTNum],dtptg[MaxGTNum];

  int ach[128],cch[128];
  float aval[128],cval[128];
  int nhit,chit;

  Float_t tthe[2],tphi[2];

  int anodemap[PPAC_NUM][PPAC_NUM];
  int cathodemap[4*PPAC_NUM][4*PPAC_NUM];
  float offTheta[PPAC_NUM],gainTheta[PPAC_NUM],quadTheta[PPAC_NUM];
  float offPhi[PPAC_NUM],gainPhi[PPAC_NUM];
  float offGT[QNum*4],gainGT[QNum*4],offdt[QNum*4];
  float offGT2[QNum*4],gainGT2[QNum*4];

  CHICOEVENT ChicoEvent;

  #ifdef SKIMCHICO
  FILE *fskim;
  SKIM Skim;
  CHICO    Chico[MaxChicoNum];
  MODE2    Mode2[MaxGTNum];
  TCutG  *etCut;
  #endif
};
