#include "StRoot/StPicoAnaTreeMaker/StAnaTree.h"
#include "StRoot/StPicoAnaTreeMaker/StPicoAnaTreeMaker.h"
#include "StRoot/StPicoAnaTreeMaker/StEventHeader.h"
#include "StRoot/StPicoAnaTreeMaker/StElectronTrack.h"
#include "StRoot/StPicoAnaTreeMaker/StMuonTrack.h"
#include "StRoot/StPicoAnaTreeMaker/StEEPair.h"
#include "StRoot/StPicoAnaTreeMaker/StEMuPair.h"
#include "StRoot/StPicoAnaTreeMaker/StMuMuPair.h"
#include "StRoot/StPicoAnaTreeMaker/StEmcTrigger.h"
#include "StRoot/StPicoAnaTreeMaker/StHadronTrack.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "mBadRunList.h"
#include "StMyAnaTreeMaker.h"

#define nCentrals 10 
#define eMass 0.000510999
#define muMass 0.1056583715

ClassImp(StMyAnaTreeMaker)

  //-----------------------------------------------------------------------------
  StMyAnaTreeMaker::StMyAnaTreeMaker(const char* name, StPicoAnaTreeMaker *treeMaker, const char* outName)
: StMaker(name)
{
  mPicoAnaTreeMaker = treeMaker;
  mAnaTree = 0;
  TH1F:: SetDefaultSumw2();//zaochen add
  mOutName = outName;

  mTrigSelect = -1;
  mHTth = 0;
  mHTAdc0th = 250;
  mEmcPtth = 2.5;

  mVzCut[0] = -100; mVzCut[1] = 100;
  mVzDiffCut[0] = -3; mVzDiffCut[1] = 3;
  mnHitsFitCut[0] = 20; mnHitsFitCut[1] = 50;
  mnHitsDedxCut[0] = 15; mnHitsDedxCut[1] = 50;
  mRatioCut[0] = 0.52; mRatioCut[1] = 1.02;
  mHFTTrackCut = false;

  mEPtCut[0] = 0.2; mEPtCut[1] = 30;
  mEEtaCut[0] = -1.; mEEtaCut[1] = 1.;
  mEDcaCut[0] = 0.; mEDcaCut[1] = 3.;
  mEInvBetaCut[0] = 0.975; mEInvBetaCut[1] = 1.025;
  mELocalYCut[0] = -2.0; mELocalYCut[1] = 1.8;
  mELocalZCut[0] = -3.05; mELocalZCut[1] = 3.05;
  mEnSigECut[0] = -0.8; mEnSigECut[1] = 2.5;
  mHTEnSigECut[0] = -1.5; mHTEnSigECut[1] = 3.;

  mEmcEPtCut[0] = 1.2; mEmcEPtCut[1] = 100;
  mEmcEEtaCut[0] = -1.; mEmcEEtaCut[1] = 1.;
  mEmcEPveCut[0] = 0.3; mEmcEPveCut[1] = 1.5;
  mEmcEDcaCut[0] = 0.; mEmcEDcaCut[1] = 3.;

  mEnEtaCut[0] = 1; mEnEtaCut[1] = 20;
  mEnPhiCut[0] = 1; mEnPhiCut[1] = 20;
  mEZDistCut[0] = -2; mEZDistCut[1] = 2;
  mEPhiDistCut[0] = -0.01; mEPhiDistCut[1] = 0.01;

  mMuPtCut[0] = 1.0; mMuPtCut[1] = 100.;
  mMuEtaCut[0] = -0.65; mMuEtaCut[1] = 0.65;
  mMuDcaCut[0] = 0.; mMuDcaCut[1] = 3.;
  mMunSigPiCut[0] = -1; mMunSigPiCut[1] = 3.;
  mMudTCut[0] = -0.7; mMudTCut[1] = 0.7;
  mMudYCut[0] = -60; mMudYCut[1] = 60;
  mMudZCut[0] = -30; mMudZCut[1] = 30;

  mDauEPtCut[0] = 0.2; mDauEPtCut[1]  = 100.;
  mDauEEtaCut[0] = -1.; mDauEEtaCut[1] = 1.;
  mDauEDcaToVtxCut[0] = 0; mDauEDcaToVtxCut[1] = 3;
  mDauEDcaDistCut[0] = 0; mDauEDcaDistCut[1] = 1;

  mDauMuPtCut[0] = 1.0; mDauMuPtCut[1] = 100.;
  mDauMuEtaCut[0] = -0.65; mDauMuEtaCut[1] = 0.65;
  mDauMuDcaToVtxCut[0] = 0.; mDauMuDcaToVtxCut[1] = 5.;
  mCosThetaStarCut[0] = -1.; mCosThetaStarCut[1] = 1.;
  mPointingAngleCut[0] = -3.14159; mPointingAngleCut[1] = 3.14159;
  mPairDcaCut[0] = 0.; mPairDcaCut[1] = 1000.;
  mPairDecayLCut[0] = 0.; mPairDecayLCut[1] = 1000.;
  mEEYCut[0] = -1.; mEEYCut[1] = 1.;
  mPairMassCut[0] = 0.; mPairMassCut[1] = 20;

  mEEYCut[0] = -1.; mEEYCut[1] = 1.;
  mEMuYCut[0] = -1.; mEMuYCut[1] = 1.;
  mMuMuYCut[0] = -0.65; mMuMuYCut[1] = 0.65;

  mPEMassCut[0] = 0.; mPEMassCut[1] = 0.02;

  mHadPtCut[0]  = 0.3;  mHadPtCut[1] = 20.0;
  mHadEtaCut[0] = -0.7; mHadEtaCut[1]= 0.7;
  mHadDcaCut[0] = 0.; mHadDcaCut[1] = 3.0;

  mNBadRuns = sizeof(mBadRuns)/sizeof(int);

  nMassBins = 500;
  massMin = 0.;
  massMax = 5.;
  nPtBins = 200;
  ptMin = 0.;
  ptMax = 10.;

  iran = 231;
  current_centrality = 0;

  fPhiVm = new TF1("phiVm","[0]/(exp(-[1]/x)+[2])",0.,0.08);
  fPhiVm->SetParameters(0.000228498,0.363774,0.000272462);

}

//----------------------------------------------------------------------------- 
StMyAnaTreeMaker::~StMyAnaTreeMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnaTreeMaker::Init() {

  if(mOutName!="") {
    fout = new TFile(mOutName.Data(),"RECREATE");
  }else{
    fout = new TFile("test.ana.root","RECREATE");
  }

  fEDcaCut = new TF1("fEDcaCut","[0]*pow(1/x,[1])+[2]",0,100);
  fEDcaCut->SetParameters(0.0121809,0.52427,0.00586945);

  fMuDcaCut = new TF1("fMuDcaCut","[0]*pow(1/x,[1])+[2]",1,100);
  fMuDcaCut->SetParameters(0.00430293,2.45282,0.00650272);

  fMudTCutLow = new TF1("fMudTCutLow","[0]*pow(1./x,[1])+[2]",1,100);
  fMudTCutLow->SetParameters(-3.62155,5.17877,-0.398502);

  fMudTCutHigh = new TF1("fMudTCutHigh","[0]*pow(1./x,[1])+[2]",1,100);
  fMudTCutHigh->SetParameters(2.5084,2.67417,0.399722);

  memset(nEventsInBuffer,0,sizeof(nEventsInBuffer));
  memset(buffer_fullFlag,0,sizeof(buffer_fullFlag));
  memset(buffer_nePlus,0,sizeof(buffer_nePlus));
  memset(buffer_neMinus,0,sizeof(buffer_neMinus));
  memset(buffer_nmuPlus,0,sizeof(buffer_nmuPlus));
  memset(buffer_nmuMinus,0,sizeof(buffer_nmuMinus));

  Int_t dblSize = sizeof(Double_t);
  memset(mMtdT0Corr,      0, 30*5*dblSize);

  ifstream inData;
  inData.open("StRoot/StMyAnaTreeMaker/Run14_AuAu200_CalibDtof.offline.dat");
  if(!inData.is_open())
  {
    LOG_ERROR << "Unable to get the T0 offset parameters from local file" <<endm;
    LOG_ERROR << "Check if this file exists: StRoot/StMyAnaTreeMaker/Run14_AuAu200_CalibDtof.offline.dat" << endm;
    return kStErr;
  }

  for(Int_t i=0;i<30;i++)
  {
    for(Int_t j=0;j<5;j++)
    {
      int backlegId, moduleId, t0Corr;
      inData >> backlegId >> moduleId >> t0Corr;
      mMtdT0Corr[backlegId-1][moduleId-1]=t0Corr;
    }
  }
  inData.close();

  declareHistograms();

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnaTreeMaker::Finish() {
  fout->cd();
  fout->Write();
  fout->Close();
  printCuts();
  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnaTreeMaker::declareHistograms() {

  fout->cd();
  hnEvents = new TH1F("hnEvents","hnEvents",10,0,10);
  hnTracks = new TH1F("hnTracks","hnTracks",10,0,10);
  hVzVpdVz = new TH2F("hVzVpdVz","hVzVpdVz;Vz (cm); Vz_{VPD} (cm);",400,-100,100,400,-100,100);
  hVzdVz = new TH2F("hVzdVz","hVzdVz;Vz (cm); Vz_{VPD}-Vz_{TPC} (cm);",400,-100,100,400,-10,10);
  hRefMultCut = new TH1F("hRefMultCut","Reference Multiplicity after cut;uncorrected dN_{ch}/d#eta;Counts",1000,0,1000);
  hVertexZCut = new TH1F("hVertexZCut","vertexZ after cut;Vz (cm);Couts",400,-100,100);

  hNe = new TH2F("hNe","#e+ vs. #e-;#e^{+} candidate;#e^{-} candidate;Counts",100,0,100,100,0,100);
  hNemu = new TH2F("hNemu","#e vs. #mu;#e candidate;#mu candidate;Counts",100,0,100,100,0,100);
  hNmu = new TH2F("hNmu","#mu+ vs. #mu-;#mu^{+} candidate;#mu^{-} candidate;Counts",100,0,100,100,0,100);

  hEEtavsPhi = new TH2F("hEEtavsPhi","hEEtavsPhi; #phi; #eta;",300,-2.0,4.5,200,-1,1);
  hEPhivsPt = new TH2F("hEPhivsPt","hEPhivsPt; q*p_{T} (GeV/c); #phi;",2.*nPtBins,-ptMax,ptMax,300,-2.0,4.5);
  hEEtavsPt = new TH2F("hEEtavsPt","hEEtavsPt; q*p_{T} (GeV/c); #eta;",2.*nPtBins,-ptMax,ptMax,200,-1,1);
  hEDcavsPt = new TH2F("hEDcavsPt","hEDcavsPt; q*p_{T} (GeV/c); dca (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);

  hEEtavsPhiwHft = new TH2F("hEEtavsPhiwHft","hEEtavsPhiwHft; #phi; #eta;",300,-2.0,4.5,200,-1,1);
  hEPhivsPtwHft = new TH2F("hEPhivsPtwHft","hEPhivsPtwHft; q*p_{T} (GeV/c); #phi;",2.*nPtBins,-ptMax,ptMax,300,-2.0,4.5);
  hEEtavsPtwHft = new TH2F("hEEtavsPtwHft","hEEtavsPtwHft; q*p_{T} (GeV/c); #eta;",2.*nPtBins,-ptMax,ptMax,200,-1,1);
  hEDcavsPtwHft = new TH2F("hEDcavsPtwHft","hEDcavsPtwHft; q*p_{T} (GeV/c); dca (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hEDcaXYvsPtwHft = new TH2F("hEDcaXYvsPtwHft","hEDcaXYvsPtwHft; q*p_{T} (GeV/c); dcaXY (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hEDcaZvsPtwHft = new TH2F("hEDcaZvsPtwHft","hEDcaZvsPtwHft; q*p_{T} (GeV/c); dcaZ (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);

  hPEUSOyOx = new TH2F("hPEUSOyOx","hPEUSOyOx; Ox (cm); Oy (cm)",600,-30,30,300,-30,30);
  hPEUSOxOz = new TH2F("hPEUSOxOz","hPEUSOxOz; Oz (cm); Ox (cm)",600,-30,30,300,-30,30);
  hPEUSOrOz = new TH2F("hPEUSOrOz","hPEUSOrOz; Oz (cm); Or (cm)",600,-30,30,300,-30,30);
  hPELSOyOx = new TH2F("hPELSOyOx","hPELSOyOx; Ox (cm); Oy (cm)",600,-30,30,300,-30,30);
  hPELSOxOz = new TH2F("hPELSOxOz","hPELSOxOz; Oz (cm); Ox (cm)",600,-30,30,300,-30,30);
  hPELSOrOz = new TH2F("hPELSOrOz","hPELSOrOz; Oz (cm); Or (cm)",600,-30,30,300,-30,30);

  hPEUSOyOxwHft = new TH2F("hPEUSOyOxwHft","hPEUSOyOxwHft; Ox (cm); Oy (cm)",600,-30,30,300,-30,30);
  hPEUSOxOzwHft = new TH2F("hPEUSOxOzwHft","hPEUSOxOzwHft; Oz (cm); Ox (cm)",600,-30,30,300,-30,30);
  hPEUSOrOzwHft = new TH2F("hPEUSOrOzwHft","hPEUSOrOzwHft; Oz (cm); Or (cm)",600,-30,30,300,-30,30);
  hPELSOyOxwHft = new TH2F("hPELSOyOxwHft","hPELSOyOxwHft; Ox (cm); Oy (cm)",600,-30,30,300,-30,30);
  hPELSOxOzwHft = new TH2F("hPELSOxOzwHft","hPELSOxOzwHft; Oz (cm); Ox (cm)",600,-30,30,300,-30,30);
  hPELSOrOzwHft = new TH2F("hPELSOrOzwHft","hPELSOrOzwHft; Oz (cm); Or (cm)",600,-30,30,300,-30,30);

  hPEUSDcavsPtwHft = new TH2F("hPEUSDcavsPtwHft","hPEUSDcavsPtwHft; q*p_{T} (GeV/c); dca (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hPEUSDcaXYvsPtwHft = new TH2F("hPEUSDcaXYvsPtwHft","hPEUSDcaXYvsPtwHft; q*p_{T} (GeV/c); dcaXY (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hPEUSDcaZvsPtwHft = new TH2F("hPEUSDcaZvsPtwHft","hPEUSDcaZvsPtwHft; q*p_{T} (GeV/c); dcaZ (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hPELSDcavsPtwHft = new TH2F("hPELSDcavsPtwHft","hPELSDcavsPtwHft; q*p_{T} (GeV/c); dca (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hPELSDcaXYvsPtwHft = new TH2F("hPELSDcaXYvsPtwHft","hPELSDcaXYvsPtwHft; q*p_{T} (GeV/c); dcaXY (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hPELSDcaZvsPtwHft = new TH2F("hPELSDcaZvsPtwHft","hPELSDcaZvsPtwHft; q*p_{T} (GeV/c); dcaZ (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);

  hPEEvPvsPt = new TH2F("hPEEvPvsPt","hPEEvPvsPt; p_{T} (GeV/c); E/p;",nPtBins,0,ptMax,500,0,5);
  hPEPvEvsPt = new TH2F("hPEPvEvsPt","hPEPvEvsPt; p_{T} (GeV/c); p/E;",nPtBins,0,ptMax,500,0,5);

  hEvPvsPt = new TH2F("hEvPvsPt","hEvPvsPt; p_{T} (GeV/c); E/p;",nPtBins,0,ptMax,500,0,5);
  hPvEvsPt = new TH2F("hPvEvsPt","hPvEvsPt; p_{T} (GeV/c); p/E;",nPtBins,0,ptMax,500,0,5);

  hEmcnSigEvsPt = new TH2F("hEmcnSigEvsPt","hEmcnSigEvsPt; p_{T} (GeV/c); n#sigma_{e};",nPtBins,0,ptMax,500,-5,5);
  hTofnSigEvsPt = new TH2F("hTofnSigEvsPt","hTofnSigEvsPt; p_{T} (GeV/c); n#sigma_{e};",nPtBins,0,ptMax,500,-5,5);

  hMuEtavsPhi = new TH2F("hMuEtavsPhi","hMuEtavsPhi; #phi; #eta;",300,-2.0,4.5,200,-1,1);
  hMuPhivsPt = new TH2F("hMuPhivsPt","hMuPhivsPt; q*p_{T} (GeV/c); #phi;",2.*nPtBins,-ptMax,ptMax,300,-2.0,4.5);
  hMuEtavsPt = new TH2F("hMuEtavsPt","hMuEtavsPt; q*p_{T} (GeV/c); #eta;",2.*nPtBins,-ptMax,ptMax,200,-1,1);
  hMuDcavsPt = new TH2F("hMuDcavsPt","hMuDcavsPt; q*p_{T} (GeV/c); dca (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);

  hMuEtavsPhiwHft = new TH2F("hMuEtavsPhiwHft","hMuEtavsPhiwHft; #phi; #eta;",300,-2.0,4.5,200,-1,1);
  hMuPhivsPtwHft = new TH2F("hMuPhivsPtwHft","hMuPhivsPtwHft; q*p_{T} (GeV/c); #phi;",2.*nPtBins,-ptMax,ptMax,300,-2.0,4.5);
  hMuEtavsPtwHft = new TH2F("hMuEtavsPtwHft","hMuEtavsPtwHft; q*p_{T} (GeV/c); #eta;",2.*nPtBins,-ptMax,ptMax,200,-1,1);
  hMuDcavsPtwHft = new TH2F("hMuDcavsPtwHft","hMuDcavsPtwHft; q*p_{T} (GeV/c); dca (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hMuDcaXYvsPtwHft = new TH2F("hMuDcaXYvsPtwHft","hMuDcaXYvsPtwHft; q*p_{T} (GeV/c); dcaXY (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);
  hMuDcaZvsPtwHft = new TH2F("hMuDcaZvsPtwHft","hMuDcaZvsPtwHft; q*p_{T} (GeV/c); dcaZ (cm);",2.*nPtBins,-ptMax,ptMax,1000,0,1);


  //hEENumInvMassvsPtMB = new TH2F("hEENumInvMassvsPtMB","same event mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  //hEEDenInvMassvsPtLikePosMB = new TH2F("hEEDenInvMassvsPtLikePosMB","like-sign pos mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  //hEEDenInvMassvsPtLikeNegMB = new TH2F("hEEDenInvMassvsPtLikeNegMB","like-sign neg mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEENumInvMassvsPtMB = new TH2F("hEENumInvMassvsPtMB","same event mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,1200,0,12);
  hEEDenInvMassvsPtLikePosMB = new TH2F("hEEDenInvMassvsPtLikePosMB","like-sign pos mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,1200,0,12);
  hEEDenInvMassvsPtLikeNegMB = new TH2F("hEEDenInvMassvsPtLikeNegMB","like-sign neg mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,1200,0,12);

  hEMuNumInvMassvsPtMB = new TH2F("hEMuNumInvMassvsPtMB","same event mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEMuDenInvMassvsPtLikePosMB = new TH2F("hEMuDenInvMassvsPtLikePosMB","like-sign pos mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEMuDenInvMassvsPtLikeNegMB = new TH2F("hEMuDenInvMassvsPtLikeNegMB","like-sign neg mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);

  hMuMuNumInvMassvsPtMB = new TH2F("hMuMuNumInvMassvsPtMB","same event mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hMuMuDenInvMassvsPtLikePosMB = new TH2F("hMuMuDenInvMassvsPtLikePosMB","like-sign pos mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hMuMuDenInvMassvsPtLikeNegMB = new TH2F("hMuMuDenInvMassvsPtLikeNegMB","like-sign neg mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);

  hEENumInvMassvsPtMBwHft = new TH2F("hEENumInvMassvsPtMBwHft","same event mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEEDenInvMassvsPtLikePosMBwHft = new TH2F("hEEDenInvMassvsPtLikePosMBwHft","like-sign pos mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEEDenInvMassvsPtLikeNegMBwHft = new TH2F("hEEDenInvMassvsPtLikeNegMBwHft","like-sign neg mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);

  hEMuNumInvMassvsPtMBwHft = new TH2F("hEMuNumInvMassvsPtMBwHft","same event mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEMuDenInvMassvsPtLikePosMBwHft = new TH2F("hEMuDenInvMassvsPtLikePosMBwHft","like-sign pos mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEMuDenInvMassvsPtLikeNegMBwHft = new TH2F("hEMuDenInvMassvsPtLikeNegMBwHft","like-sign neg mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);

  hMuMuNumInvMassvsPtMBwHft = new TH2F("hMuMuNumInvMassvsPtMBwHft","same event mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hMuMuDenInvMassvsPtLikePosMBwHft = new TH2F("hMuMuDenInvMassvsPtLikePosMBwHft","like-sign pos mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hMuMuDenInvMassvsPtLikeNegMBwHft = new TH2F("hMuMuDenInvMassvsPtLikeNegMBwHft","like-sign neg mass spectrum in MBwHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);


  hEENumInvMassvsPtMBnophiv = new TH2F("hEENumInvMassvsPtMBnophiv","same event mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEEDenInvMassvsPtLikePosMBnophiv = new TH2F("hEEDenInvMassvsPtLikePosMBnophiv","like-sign pos mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
  hEEDenInvMassvsPtLikeNegMBnophiv = new TH2F("hEEDenInvMassvsPtLikeNegMBnophiv","like-sign neg mass spectrum in MB;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);


  hUSphivM = new TH2F("hUSphivM","unlike-sign phiv vs m_{ee} ;M_{ee} (GeV/c^{2});#phi_{V} ",3200,0,3.2,1000,0,3.2);
  hLSPosphivM = new TH2F("hLSPosphivM","like-sign-pos phiv vs m_{ee} ;M_{ee} (GeV/c^{2});#phi_{V} ",3200,0,3.2,1000,0,3.2);
  hLSNegphivM = new TH2F("hLSNegphivM","like-sign-neg phiv vs m_{ee} ;M_{ee} (GeV/c^{2});#phi_{V} ",3200,0,3.2,1000,0,3.2);


  char name[200],title[200];
  sprintf(name,"hEENumInvMassvsPt");
  sprintf(title,"same event mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hEENumInvMassvsPt = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  sprintf(name,"hEEDenInvMassvsPtLikePos");
  sprintf(title,"like-sign pos mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hEEDenInvMassvsPtLikePos = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  sprintf(name,"hEEDenInvMassvsPtLikeNeg");
  sprintf(title,"like-sign neg mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hEEDenInvMassvsPtLikeNeg = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  sprintf(name,"hEMuNumInvMassvsPt");
  sprintf(title,"same event mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hEMuNumInvMassvsPt = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  sprintf(name,"hEMuDenInvMassvsPtLikePos");
  sprintf(title,"like-sign pos mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hEMuDenInvMassvsPtLikePos = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  sprintf(name,"hEMuDenInvMassvsPtLikeNeg");
  sprintf(title,"like-sign neg mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hEMuDenInvMassvsPtLikeNeg = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);


  sprintf(name,"hMuMuNumInvMassvsPt");
  sprintf(title,"same event mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hMuMuNumInvMassvsPt = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  sprintf(name,"hMuMuDenInvMassvsPtLikePos");
  sprintf(title,"like-sign pos mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hMuMuDenInvMassvsPtLikePos = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  sprintf(name,"hMuMuDenInvMassvsPtLikeNeg");
  sprintf(title,"like-sign neg mass spectrum;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality");
  hMuMuDenInvMassvsPtLikeNeg = new TH3F(name, title,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hEENumInvMassvsPtwHft = new TH3F("hEENumInvMassvsPtwHft","same event mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2});Centrality",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEEDenInvMassvsPtLikePoswHft = new TH3F("hEEDenInvMassvsPtLikePoswHft","like-sign pos mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEEDenInvMassvsPtLikeNegwHft = new TH3F("hEEDenInvMassvsPtLikeNegwHft","like-sign neg mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hEMuNumInvMassvsPtwHft = new TH3F("hEMuNumInvMassvsPtwHft","same event mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEMuDenInvMassvsPtLikePoswHft = new TH3F("hEMuDenInvMassvsPtLikePoswHft","like-sign pos mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEMuDenInvMassvsPtLikeNegwHft = new TH3F("hEMuDenInvMassvsPtLikeNegwHft","like-sign neg mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hMuMuNumInvMassvsPtwHft = new TH3F("hMuMuNumInvMassvsPtwHft","same event mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hMuMuDenInvMassvsPtLikePoswHft = new TH3F("hMuMuDenInvMassvsPtLikePoswHft","like-sign pos mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hMuMuDenInvMassvsPtLikeNegwHft = new TH3F("hMuMuDenInvMassvsPtLikeNegwHft","like-sign neg mass spectrum in wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hEEDenInvMassvsPtMix = new TH3F("hEEDenInvMassvsPtMix","mixed unlike-sign event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEEDenInvMassvsPtMixLikePos = new TH3F("hEEDenInvMassvsPtMixLikePos","mixed like-sign ++ event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEEDenInvMassvsPtMixLikeNeg = new TH3F("hEEDenInvMassvsPtMixLikeNeg","mixed like-sign -- event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hEMuDenInvMassvsPtMix = new TH3F("hEMuDenInvMassvsPtMix","mixed unlike-sign event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEMuDenInvMassvsPtMixLikePos = new TH3F("hEMuDenInvMassvsPtMixLikePos","mixed like-sign ++ event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEMuDenInvMassvsPtMixLikeNeg = new TH3F("hEMuDenInvMassvsPtMixLikeNeg","mixed like-sign -- event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hMuMuDenInvMassvsPtMix = new TH3F("hMuMuDenInvMassvsPtMix","mixed unlike-sign event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hMuMuDenInvMassvsPtMixLikePos = new TH3F("hMuMuDenInvMassvsPtMixLikePos","mixed like-sign ++ event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hMuMuDenInvMassvsPtMixLikeNeg = new TH3F("hMuMuDenInvMassvsPtMixLikeNeg","mixed like-sign -- event mass spectrum ;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  //wHFT
  hEEDenInvMassvsPtMixwHft = new TH3F("hEEDenInvMassvsPtMixwHft","mixed unlike-sign event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEEDenInvMassvsPtMixLikePoswHft = new TH3F("hEEDenInvMassvsPtMixLikePoswHft","mixed like-sign ++ event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEEDenInvMassvsPtMixLikeNegwHft = new TH3F("hEEDenInvMassvsPtMixLikeNegwHft","mixed like-sign -- event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hEMuDenInvMassvsPtMixwHft = new TH3F("hEMuDenInvMassvsPtMixwHft","mixed unlike-sign event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEMuDenInvMassvsPtMixLikePoswHft = new TH3F("hEMuDenInvMassvsPtMixLikePoswHft","mixed like-sign ++ event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hEMuDenInvMassvsPtMixLikeNegwHft = new TH3F("hEMuDenInvMassvsPtMixLikeNegwHft","mixed like-sign -- event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  hMuMuDenInvMassvsPtMixwHft = new TH3F("hMuMuDenInvMassvsPtMixwHft","mixed unlike-sign event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hMuMuDenInvMassvsPtMixLikePoswHft = new TH3F("hMuMuDenInvMassvsPtMixLikePoswHft","mixed like-sign ++ event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);
  hMuMuDenInvMassvsPtMixLikeNegwHft = new TH3F("hMuMuDenInvMassvsPtMixLikeNegwHft","mixed like-sign -- event mass spectrum wHft;p_{T} (GeV/c);M_{inv} (GeV/c^{2})",nPtBins,ptMin,ptMax,nMassBins,massMin,massMax,nCentrals,0,nCentrals);

  cout << "Declare Had Hists" << endl;
  hHadPt  = new TH1F("hHadPt","Hadron Branch Pt Spectrum; p_{T} (GeV/c); Counts",nPtBins,ptMin,ptMax);
  hEPt  = new TH1F("hEPt","Electron Branch Pt Spectrum; p_{T} (GeV/c); Counts",nPtBins,ptMin,ptMax);
  hEEPt_US  = new TH1F("hEEPt_US","Electron Pair US Pt Spectrum; p_{T} (GeV/c); Counts",nPtBins,ptMin,ptMax);
  hEEPt_LS  = new TH1F("hEEPt_LS","Electron Pair LS Pt Spectrum; p_{T} (GeV/c); Counts",nPtBins,ptMin,ptMax);
  hHadDca = new TH1F("hHadDca","Hadron Branch DCA; DCA (cm); Counts",1000,0,1); 
  hHadDcaXY = new TH1F("hHadDcaXY","Hadron Branch DCAXY; DCA_XY (cm); Counts",1000,0,1); 
  hHadDcaZ = new TH1F("hHadDcaZ","Hadron Branch DCAZ; DCA_Z (cm); Counts",1000,0,1); 
  hHadPhi = new TH1F("hHadPhi","Hadron Branch Phi; Phi; Counts",300,-3.2,3.2); 

  hHadEtaPhi = new TH2F("hHadEtaPhi","Hadron Branch Phi vs Pseudorapidity; Eta; Phi",200,-1,1,300,-3.2,3.2); 
  hHadEDelPhiPt = new TH2F("hHadEDelPhiPt","Hadron-Electron DeltaPhi vs Pt; p_{T} (GeV/c); #Delta#phi (rads)",nPtBins,ptMin,ptMax,300,-2.0,4.5);
  hHadEDelPhiPt_high = new TH2F("hHadEDelPhiPt_high","Hadron-Electron DeltaPhi vs Electron Pt (hp_{T} > 2 GeV/c^{2}); p_{T} (GeV/c); #Delta#phi (rads)",nPtBins,ptMin,ptMax,300,-2.0,4.5);
  hHadEEDelPhiPt_LS = new TH2F("hHadEEDelPhiPt_LS","Hadron-EE LS DeltaPhi vs Pt; p_{T} (GeV/c); #Delta#phi (rads)",nPtBins,ptMin,ptMax,300,-2.0,4.5);
  hHadEEDelPhiPt_US = new TH2F("hHadEEDelPhiPt_US","Hadron-EE US DeltaPhi vs Pt; p_{T} (GeV/c); #Delta#phi (rads)",nPtBins,ptMin,ptMax,300,-2.0,4.5);
  hHadMuDelPhiPt = new TH2F("hHadMuDelPhiPt","Hadron-Muon DeltaPhi vs Pt; p_{T} (GeV/c); #Delta#phi (rads)",nPtBins,ptMin,ptMax,300,-2.0,4.5);
  hHadDcaPt = new TH2F("hHadDcaPt","Hadron Branch DCA vs pT; q*p_{T} (GeV/c); DCA (cm)",2*nPtBins,-ptMax,ptMax,1000,0,3); 
  hHadPtEPt = new TH2F("hHadPtEPt","Hadron pT vs Electron pT; Electron p_{T} (GeV/c); Hadron p_{T} (GeV/c)",nPtBins,ptMin,ptMax,nPtBins,ptMin,ptMax); 
}


//----------------------------------------------------------------------------- 
void StMyAnaTreeMaker::Clear(Option_t *opt) {

}

//----------------------------------------------------------------------------- 
Int_t StMyAnaTreeMaker::Make() {

  if(!mPicoAnaTreeMaker) {
    LOG_WARN << " No PicoAnaTreeMaker! Skip! " << endm;
    return kStOK;
  }

  mAnaTree = mPicoAnaTreeMaker->anaTree();

  if(!mAnaTree) {
    LOG_WARN << " No AnaTree! Skip! " << endm;
    return kStWarn;
  }

  hnEvents->Fill(0);
  int runId = mAnaTree->event()->runId();
  for(int i=0;i<mNBadRuns;i++){
    if(runId==mBadRuns[i]) return kStOK;
  }
  hnEvents->Fill(1);
  int eventPass = 0;
  if(mTrigSelect<0) eventPass = 1;
  else if(mTrigSelect==0&&mAnaTree->event()->isMinBias()) eventPass = 1;
  else if(mTrigSelect==1&&mAnaTree->event()->isHT11()){ eventPass = 1; mHTth = 11; mHTAdc0th = 180; mEmcPtth = 1.5;} //HT0 180
  else if(mTrigSelect==2&&mAnaTree->event()->isHT15()){ eventPass = 1; mHTth = 15; mHTAdc0th = 250; mEmcPtth = 2.5;} //HT1 240
  else if(mTrigSelect==3&&mAnaTree->event()->isHT18()){ eventPass = 1;mHTth = 18; mHTAdc0th = 300; mEmcPtth = 3.0;}  //HT2 300
  else if(mTrigSelect==4&&mAnaTree->event()->isHT25()){ eventPass = 1;mHTth = 25; mHTAdc0th = 400; mEmcPtth = 4.0;}  //HT3 400
  else if(mTrigSelect==5&&mAnaTree->event()->isEMuon()){ eventPass = 1;mHTth = 13; mHTAdc0th = 210; mEmcPtth = 2.;} //EMu 210
  else if(mTrigSelect==6&&mAnaTree->event()->isDiMuon()) eventPass = 1; //di-muon
  else if(mTrigSelect==7&&mAnaTree->event()->isSingleMuon()) eventPass = 1; //single-muon

  if(eventPass==0) return kStOK;
  hnEvents->Fill(2);

  Double_t vzVpd=mAnaTree->event()->vzVpd();
  StThreeVectorF mPrimaryVertex = mAnaTree->event()->primaryVertex();
  hVzVpdVz->Fill(mPrimaryVertex.z(),vzVpd);
  hVzdVz->Fill(mPrimaryVertex.z(), vzVpd-mPrimaryVertex.z());
  if(mPrimaryVertex.z()<mVzCut[0]||mPrimaryVertex.z()>mVzCut[1]) return kStOK;
  hnEvents->Fill(3);
  if(vzVpd-mPrimaryVertex.z()<mVzDiffCut[0]||vzVpd-mPrimaryVertex.z()>mVzDiffCut[1]) return kStOK;
  hnEvents->Fill(4);

  hRefMultCut->Fill(mAnaTree->event()->grefMult());
  hVertexZCut->Fill(mPrimaryVertex.z());

  int centrality = getCentrality();
  current_centrality = getCentrality();

  int magBufferPointer = -1;
  double bfield = mAnaTree->event()->bField();
  if(bfield>0) magBufferPointer = 0; //ff
  if(bfield<0) magBufferPointer = 1;//rff 
  if(magBufferPointer<0) return kStOK;
  int cenBufferPointer = centrality-1;
  if(cenBufferPointer<0||cenBufferPointer>=nCenBin) return kStOK;
  double current_eventplane = mAnaTree->event()->eventPlane();
  int   eveBufferPointer = int((current_eventplane-0.)/TMath::Pi()*nEveBin);
  if(eveBufferPointer<0 ||eveBufferPointer>=nEveBin) return kStOK;
  double Vz = mAnaTree->event()->primaryVertex().z();
  int vzBufferPointer=int((Vz-mVzCut[0])/(mVzCut[1]-mVzCut[0])*nVzBin);


  int current_nePlus = 0, current_neMinus = 0;
  int current_nmuPlus = 0, current_nmuMinus = 0;
  int nE = mAnaTree->numberOfETracks();

  memset(current_ePlusFlag,0,sizeof(current_ePlusFlag));
  memset(current_eMinusFlag,0,sizeof(current_eMinusFlag));
  memset(current_ePlusIsHFT,0,sizeof(current_ePlusIsHFT));
  memset(current_eMinusIsHFT,0,sizeof(current_eMinusIsHFT));
  for(int i=0;i<nE;i++){
    StElectronTrack *eTrk = (StElectronTrack*)mAnaTree->eTrack(i);
    int charge = eTrk->charge();
    int eHTflag = 0, eflag = 0, etrgflag = 0;
    if(passHTEIDCuts(eTrk)) eHTflag = 1;
    if(passEIDCuts(eTrk)) eflag = 1;
    if(isHTTrigE(eTrk)&&eHTflag) etrgflag = 1;

    if(!eHTflag&&!eflag) continue;
    double pt = eTrk->pMom().perp();
    double eta = eTrk->pMom().pseudoRapidity();
    double phi = eTrk->pMom().phi();
    if(eHTflag){
      if(eTrk->pve()>0) hEvPvsPt->Fill(pt,1./eTrk->pve());
      hPvEvsPt->Fill(pt,eTrk->pve());
      hEmcnSigEvsPt->Fill(pt,eTrk->nSigmaElectron());
    }
    if(eflag) hTofnSigEvsPt->Fill(pt,eTrk->nSigmaElectron());
    if(charge>0){ 
      current_ePlusFlag[current_nePlus] = eflag + 2*eHTflag+4*etrgflag;
      if(eTrk->isHFTTrack()) current_ePlusIsHFT[current_nePlus] = 1;
      current_ePlus[current_nePlus].SetPtEtaPhiM(pt,eta,phi,eMass);
      current_nePlus++; 
      hnTracks->Fill(1); 
    }
    if(charge<0){ 
      current_eMinusFlag[current_neMinus] = eflag + 2*eHTflag+4*etrgflag;
      if(eTrk->isHFTTrack()) current_eMinusIsHFT[current_neMinus] = 1;
      current_eMinus[current_neMinus].SetPtEtaPhiM(pt,eta,phi,eMass);
      current_neMinus++; 
      hnTracks->Fill(2); 
    }
    double p = eTrk->pMom().mag();
    double dca = eTrk->dca();
    double dcaXY = eTrk->dcaXY();
    double dcaZ = eTrk->dcaZ();
    int isHft = eTrk->isHFTTrack();

    hEEtavsPhi->Fill(phi,eta);
    hEPhivsPt->Fill(pt*charge,phi);
    hEEtavsPt->Fill(pt*charge,eta);
    hEDcavsPt->Fill(pt*charge,dca);
    if(isHft){
      hEEtavsPhiwHft->Fill(phi,eta);
      hEPhivsPtwHft->Fill(pt*charge,phi);
      hEEtavsPtwHft->Fill(pt*charge,eta);
      hEDcavsPtwHft->Fill(pt*charge,dca);
      hEDcaXYvsPtwHft->Fill(pt*charge,dcaXY);
      hEDcaZvsPtwHft->Fill(pt*charge,dcaZ);
    }

    // Hadron Comparisons
    hEPt->Fill(pt); // need for normalization
    int nHad = mAnaTree->numberOfHTracks(); 
    for(int j=0; j < nHad; j++)
    {
      StHadronTrack *hTrk = (StHadronTrack*)mAnaTree->hTrack(j);
      if(!passHadronCuts(hTrk)|| hTrk->id() == eTrk->id()) continue;
      int charge = hTrk->charge();
      double hpt = hTrk->gMom().perp();
      double heta = hTrk->gMom().pseudoRapidity();
      double hphi = hTrk->gMom().phi();
      double hdca = hTrk->dca();
      double hdcaXY = hTrk->dcaXY();
      double hdcaZ  = hTrk->dcaZ();

      hHadPt    ->Fill(hpt);
      hHadPhi   ->Fill(hphi);
      hHadDca   ->Fill(hdca);
      hHadDcaXY ->Fill(hdcaXY);
      hHadDcaZ  ->Fill(hdcaZ);
      hHadEtaPhi->Fill(heta,hphi);
      hHadDcaPt ->Fill(charge*hpt,hdca);
      hHadPtEPt ->Fill(pt,hpt);

      double dphi = hphi - phi;
      dphi = delPhiCorrect(dphi);
      if(pt > 1.0 && hpt > 1.5 && hdca < 1.0) hHadEDelPhiPt->Fill(pt,dphi); // electron pt
      if(pt > 1.0 && hpt > 4.0 && hdca < 1.0) hHadEDelPhiPt_high->Fill(pt,dphi); // electron pt
    }
  }

  memset(current_muPlusFlag,0,sizeof(current_muPlusFlag));
  memset(current_muMinusFlag,0,sizeof(current_muMinusFlag));
  memset(current_muPlusIsHFT,0,sizeof(current_muPlusIsHFT));
  memset(current_muMinusIsHFT,0,sizeof(current_muMinusIsHFT));
  int nMu = mAnaTree->numberOfMuTracks();
  for(int i=0;i<nMu;i++){
    StMuonTrack *muTrk = (StMuonTrack*)mAnaTree->muTrack(i);
    int charge = muTrk->charge();
    if(!passMuIDCuts(muTrk));
    double pt = muTrk->pMom().perp();
    double eta = muTrk->pMom().pseudoRapidity();
    double phi = muTrk->pMom().phi();
    if(charge>0){ 
      current_muPlus[current_nmuPlus].SetPtEtaPhiM(pt,eta,phi,eMass);
      if(muTrk->triggerFlag()>0) current_muPlusFlag[current_nmuPlus] = 1;
      if(muTrk->isHFTTrack()) current_muPlusIsHFT[current_nmuPlus] = 1;
      current_nmuPlus++; 
      hnTracks->Fill(3);
    } 
    if(charge<0){ 
      current_muMinus[current_nmuMinus].SetPtEtaPhiM(pt,eta,phi,eMass);
      if(muTrk->triggerFlag()>0) current_muMinusFlag[current_nmuPlus] = 1;
      if(muTrk->isHFTTrack()) current_muMinusIsHFT[current_nmuMinus] = 1;
      current_nmuMinus++; 
      hnTracks->Fill(4);
    }
    double p = muTrk->pMom().mag();
    double dca = muTrk->dca();
    double dcaXY = muTrk->dcaXY();
    double dcaZ = muTrk->dcaZ();
    int isHft = muTrk->isHFTTrack();

    hMuEtavsPhi->Fill(phi,eta);
    hMuPhivsPt->Fill(pt*charge,phi);
    hMuEtavsPt->Fill(pt*charge,eta);
    hMuDcavsPt->Fill(pt*charge,dca);
    if(isHft){
      hMuEtavsPhiwHft->Fill(phi,eta);
      hMuPhivsPtwHft->Fill(pt*charge,phi);
      hMuEtavsPtwHft->Fill(pt*charge,eta);
      hMuDcavsPtwHft->Fill(pt*charge,dca);
      hMuDcaXYvsPtwHft->Fill(pt*charge,dcaXY);
      hMuDcaZvsPtwHft->Fill(pt*charge,dcaZ);
    }
    // Loop over hadrons for muon correlation
    int nHad = mAnaTree->numberOfHTracks(); 
    for(int j=0; j < nHad; j++)
    {
      StHadronTrack *hTrk = (StHadronTrack*)mAnaTree->hTrack(j);
      if(!passHadronCuts(hTrk) || hTrk->id() == muTrk->id()) continue;
      int charge = hTrk->charge();
      double hpt = hTrk->gMom().perp();
      double heta = hTrk->gMom().pseudoRapidity();
      double hphi = hTrk->gMom().phi();
      double hdca = hTrk->dca();
      double hdcaXY = hTrk->dcaXY();
      double hdcaZ  = hTrk->dcaZ();

      hHadPt    ->Fill(hpt);
      hHadPhi   ->Fill(hphi);
      hHadDca   ->Fill(hdca);
      hHadDcaXY ->Fill(hdcaXY);
      hHadDcaZ  ->Fill(hdcaZ);
      hHadEtaPhi->Fill(heta,hphi);

      double dphi = hphi - phi;
      dphi = delPhiCorrect(dphi);
      hHadMuDelPhiPt->Fill(pt,dphi); // Muon pt
    }
  }


  int nEEPairs = mAnaTree->numberOfEEPairs();
  for(int iee=0; iee< nEEPairs;iee++){	
    StEEPair* ee = (StEEPair*) mAnaTree->eePair(iee);
    int dauIndex1 = ee->dauIndex1();
    int dauIndex2 = ee->dauIndex2();
    StElectronTrack *eTrk1 = mAnaTree->eTrack(dauIndex1);
    StElectronTrack *eTrk2 = mAnaTree->eTrack(dauIndex2);
    int flag = 0;
    if(mTrigSelect==0){
      if(passEIDCuts(eTrk1)&&passEIDCuts(eTrk2)) flag=1;
    }
    if(mTrigSelect==1||mTrigSelect==2||mTrigSelect==3||mTrigSelect==4||mTrigSelect==5){
      int htTrkFlag = 0,eTrkFlag = 0;
      if(isHTTrigE(eTrk1)&&passHTEIDCuts(eTrk1)) htTrkFlag += 1;
      //if(isHTTrigE(eTrk1)&&isHTTrigE(eTrk2)) continue;
      if(isHTTrigE(eTrk2)&&passHTEIDCuts(eTrk2)) htTrkFlag +=1;
      //if(passEIDCuts(eTrk1)||passHTEIDCuts(eTrk1)) eTrkFlag++;
      if(passHTEIDCuts(eTrk1)&&passHTEIDCuts(eTrk2)) eTrkFlag++;
      //if(htTrkFlag>=1&&eTrkFlag==2) flag = 1;
      if(htTrkFlag>=1&&eTrkFlag>=1) flag = 1;
    }
    if(mTrigSelect==-1){
      if(passHTEIDCuts(eTrk1)&&passHTEIDCuts(eTrk2)) flag=1;
    }
    if(flag==0) continue;

    StThreeVectorF mom1 = eTrk1->pMom();
    StThreeVectorF mom2 = eTrk2->pMom();
    StLorentzVectorF dau1(mom1,mom1.massHypothesis(eMass));
    StLorentzVectorF dau2(mom2,mom2.massHypothesis(eMass));

    double p1 = mom1.mag();
    double p2 = mom2.mag();
    double pHigh = p1>p2?p1:p2;
    double pLow = p1>p2?p2:p1;

    //check HLT
    if(pHigh<1.4||pLow<1.2) continue;

    StLorentzVectorF pair;
    pair = dau1+dau2;

    double pt = pair.perp();
    double y = pair.rapidity();
    if(!passEEPairCuts(y)) continue;
    double pmass = pair.m();
    double mass = pair.m();
    double phiV = ee->pairPhiV();
    double dauDcaDist = ee->pairDca();
    int dauIsHft1 = eTrk1->isHFTTrack();
    int dauIsHft2 = eTrk2->isHFTTrack();
    int charge1 = eTrk1->charge();
    int charge2 = eTrk2->charge();

    //if(centrality<=0) continue;
    if(charge1!=charge2){ 
      //if(pt>=0.8&&pt<0.85&&mass>=1.45&&mass<1.455){
      //	cout<<"pt = "<<pt<<" mass = "<<mass<<" y = "<<y<<endl;
      //	cout<<"index1 = "<<dauIndex1<<endl;
      //	eTrk1->Print();
      //	cout<<"index2 = "<<dauIndex2<<endl;
      //	eTrk2->Print();
      //}
      hEENumInvMassvsPtMB->Fill(pt,mass);
      hEENumInvMassvsPt->Fill(pt,mass,centrality);
      if(dauIsHft1&&dauIsHft2){ 
        hEENumInvMassvsPtMBwHft->Fill(pt,mass);
        hEENumInvMassvsPtwHft->Fill(pt,mass,centrality);
      }
    }
    if(charge1==1&&charge2==1){ 
      hEEDenInvMassvsPtLikePosMB->Fill(pt,mass);
      hEEDenInvMassvsPtLikePos->Fill(pt,mass,centrality);
      if(dauIsHft1&&dauIsHft2){ 
        hEEDenInvMassvsPtLikePosMBwHft->Fill(pt,mass);
        hEEDenInvMassvsPtLikePoswHft->Fill(pt,mass,centrality);
      }
    }
    if(charge1==-1&&charge2==-1){ 
      hEEDenInvMassvsPtLikeNegMB->Fill(pt,mass);
      hEEDenInvMassvsPtLikeNeg->Fill(pt,mass,centrality);
      if(dauIsHft1&&dauIsHft2){ 
        hEEDenInvMassvsPtLikeNegMBwHft->Fill(pt,mass);
        hEEDenInvMassvsPtLikeNegwHft->Fill(pt,mass,centrality);
      }
    }

    //photonic electron
    if(charge1!=charge2){
      hUSphivM->Fill(mass,phiV);
    }else{
      if(charge1==1)  hLSPosphivM->Fill(mass,phiV);
      if(charge1==-1) hLSNegphivM->Fill(mass,phiV);
    }
    double vcut = fPhiVm->Eval(mass);
    StThreeVectorF origin = ee->pairOrigin();
    if(pmass>mPEMassCut[0]&&pmass<mPEMassCut[1]){//&&phiV<vcut){
      //if(pmass>mPEMassCut[0]&&pmass<mPEMassCut[1]&&dauDcaDist>mDauEDcaDistCut[0]&&dauDcaDist<mDauEDcaDistCut[1])
      double pt1 = eTrk1->pMom().perp();
      double dca1 = eTrk1->dca();
      double dcaXY1 = eTrk1->dcaXY();
      double dcaZ1 = eTrk1->dcaZ();
      int    isHft1 = eTrk1->isHFTTrack();

      double pt2 = eTrk2->pMom().perp();
      double dca2 = eTrk2->dca();
      double dcaXY2 = eTrk2->dcaXY();
      double dcaZ2 = eTrk2->dcaZ();
      int    isHft2 = eTrk2->isHFTTrack();

      float pve1 = eTrk1->pve();
      float pve2 = eTrk2->pve();
      float evp1 = pve1==0?0:1./pve1;
      float evp2 = pve2==0?0:1./pve2;

      if(charge1!=charge2){
        hPEUSOyOx->Fill(origin.x(),origin.y());
        hPEUSOxOz->Fill(origin.z(),origin.x());
        hPEUSOrOz->Fill(origin.z(),origin.perp());

        if(passHTEIDCuts(eTrk1)){ 
          hPEEvPvsPt->Fill(pt1,evp1);
          hPEPvEvsPt->Fill(pt1,pve1);
        }
        if(passHTEIDCuts(eTrk2)){ 
          hPEEvPvsPt->Fill(pt2,evp2);
          hPEPvEvsPt->Fill(pt2,pve2);
        }
        if(isHft1&&isHft2){
          hPEUSOyOxwHft->Fill(origin.x(),origin.y());
          hPEUSOxOzwHft->Fill(origin.z(),origin.x());
          hPEUSOrOzwHft->Fill(origin.z(),origin.perp());
          hPEUSDcavsPtwHft->Fill(pt1*charge1,dca1);
          hPEUSDcaXYvsPtwHft->Fill(pt1*charge1,dcaXY1);
          hPEUSDcaZvsPtwHft->Fill(pt1*charge1,dcaZ1);
          hPEUSDcavsPtwHft->Fill(pt2*charge2,dca2);
          hPEUSDcaXYvsPtwHft->Fill(pt2*charge2,dcaXY2);
          hPEUSDcaZvsPtwHft->Fill(pt2*charge2,dcaZ2);
        }

      }else{
        hPELSOyOx->Fill(origin.x(),origin.y());
        hPELSOxOz->Fill(origin.z(),origin.x());
        hPELSOrOz->Fill(origin.z(),origin.perp());

        if(isHft1&&isHft2){
          hPELSOyOxwHft->Fill(origin.x(),origin.y());
          hPELSOxOzwHft->Fill(origin.z(),origin.x());
          hPELSOrOzwHft->Fill(origin.z(),origin.perp());
          hPELSDcavsPtwHft->Fill(pt1*charge1,dca1);
          hPELSDcaXYvsPtwHft->Fill(pt1*charge1,dcaXY1);
          hPELSDcaZvsPtwHft->Fill(pt1*charge1,dcaZ1);
          hPELSDcavsPtwHft->Fill(pt2*charge2,dca2);
          hPELSDcaXYvsPtwHft->Fill(pt2*charge2,dcaXY2);
          hPELSDcaZvsPtwHft->Fill(pt2*charge2,dcaZ2);
        }
      }
    }

    // Loop over hadrons for Eh correlation
    double phi1 = eTrk1->pMom().phi();
    double phi2 = eTrk2->pMom().phi();
    double pt1 = eTrk1->pMom().perp();
    double pt2 = eTrk2->pMom().perp();

    int nHad = mAnaTree->numberOfHTracks(); 
    for(int j=0; j < nHad; j++)
    {
      StHadronTrack *hTrk = (StHadronTrack*)mAnaTree->hTrack(j);
      if(!passHadronCuts(hTrk) || hTrk->id() == eTrk1->id() || hTrk->id() == eTrk2->id()) continue;
      int charge = hTrk->charge();
      double hpt = hTrk->gMom().perp();
      double heta = hTrk->gMom().pseudoRapidity();
      double hphi = hTrk->gMom().phi();
      double hdca = hTrk->dca();
      double hdcaXY = hTrk->dcaXY();
      double hdcaZ  = hTrk->dcaZ();

      double dphi = hphi - phi1;
      dphi = delPhiCorrect(dphi);
      if(charge1==charge2) 
      {
        hEEPt_LS->Fill(pt1); // need for normalization
        hHadEEDelPhiPt_LS->Fill(pt1,dphi); // Electron pt
      }
      if(charge1!=charge2) 
      {
        hEEPt_US->Fill(pt1); // need for normalization
        hHadEEDelPhiPt_US->Fill(pt1,dphi); // Muon pt
      }
    }
    //} //original end of p/E mass cut
  }

  int nEMuPairs = mAnaTree->numberOfEMuPairs();
  for(int iemu=0; iemu< nEMuPairs;iemu++){
    StEMuPair* emu = (StEMuPair*) mAnaTree->eMuPair(iemu);
    int dauIndex1 = emu->dauIndex1();
    int dauIndex2 = emu->dauIndex2();
    StElectronTrack *eTrk1 = mAnaTree->eTrack(dauIndex1);
    StMuonTrack *muTrk2 = mAnaTree->muTrack(dauIndex2);
    if(mTrigSelect==0){
      if(!passEIDCuts(eTrk1)) continue;
    }
    if(mTrigSelect==1||mTrigSelect==2||mTrigSelect==3||mTrigSelect==4||mTrigSelect==5){
      if(!passHTEIDCuts(eTrk1)) continue;
      if(!isHTTrigE(eTrk1)) continue;
    }
    if(!passMuIDCuts(muTrk2)) continue;

    StThreeVectorF mom1 = eTrk1->pMom();
    StThreeVectorF mom2 = muTrk2->pMom();
    StLorentzVectorF dau1(mom1,mom1.massHypothesis(eMass));
    StLorentzVectorF dau2(mom2,mom2.massHypothesis(muMass));

    StLorentzVectorF pair;
    pair = dau1+dau2;

    double pt = pair.perp();
    double y = pair.rapidity();
    if(!passEMuPairCuts(y)) continue;
    double mass = pair.m();
    int dauIsHft1 = eTrk1->isHFTTrack();
    int dauIsHft2 = muTrk2->isHFTTrack();
    int charge1 = eTrk1->charge();
    int charge2 = muTrk2->charge();

    double dca1 = eTrk1->dca();
    double dca2 = muTrk2->dca();

    double pt1 = eTrk1->pMom().perp();
    double pt2 = muTrk2->pMom().perp();

    //if(centrality<=0) continue;
    if(charge1!=charge2){ 
      hEMuNumInvMassvsPtMB->Fill(pt,mass);
      hEMuNumInvMassvsPt->Fill(pt,mass,centrality);
      if(dauIsHft1&&dauIsHft2&&dca1>fEDcaCut->Eval(pt1)&&dca2>fMuDcaCut->Eval(pt2)){ 
        hEMuNumInvMassvsPtMBwHft->Fill(pt,mass);
        hEMuNumInvMassvsPtwHft->Fill(pt,mass,centrality);
      }
    }
    if(charge1==1&&charge2==1){ 
      hEMuDenInvMassvsPtLikePosMB->Fill(pt,mass);
      hEMuDenInvMassvsPtLikePos->Fill(pt,mass,centrality);
      //if(dauIsHft1&&dauIsHft2) 
      if(dauIsHft1&&dauIsHft2&&dca1>fEDcaCut->Eval(pt1)&&dca2>fMuDcaCut->Eval(pt2)){ 
        hEMuDenInvMassvsPtLikePosMBwHft->Fill(pt,mass);
        hEMuDenInvMassvsPtLikePoswHft->Fill(pt,mass,centrality);
      }
    }
    if(charge1==-1&&charge2==-1){ 
      hEMuDenInvMassvsPtLikeNegMB->Fill(pt,mass);
      hEMuDenInvMassvsPtLikeNeg->Fill(pt,mass,centrality);
      //if(dauIsHft1&&dauIsHft2) 
      if(dauIsHft1&&dauIsHft2&&dca1>fEDcaCut->Eval(pt1)&&dca2>fMuDcaCut->Eval(pt2)){ 
        hEMuDenInvMassvsPtLikeNegMBwHft->Fill(pt,mass);
        hEMuDenInvMassvsPtLikeNegwHft->Fill(pt,mass,centrality);
      }
    }
  }

  int nMuMuPairs = mAnaTree->numberOfMuMuPairs();
  for(int imumu=0; imumu< nMuMuPairs;imumu++){
    StMuMuPair* mumu = (StMuMuPair*) mAnaTree->MuMuPair(imumu);
    int dauIndex1 = mumu->dauIndex1();
    int dauIndex2 = mumu->dauIndex2();
    StMuonTrack *muTrk1 = mAnaTree->muTrack(dauIndex1);
    StMuonTrack *muTrk2 = mAnaTree->muTrack(dauIndex2);
    if(!passMuIDCuts(muTrk1)) continue;
    if(!passMuIDCuts(muTrk2)) continue;

    StThreeVectorF mom1 = muTrk1->pMom();
    StThreeVectorF mom2 = muTrk2->pMom();
    StLorentzVectorF dau1(mom1,mom1.massHypothesis(muMass));
    StLorentzVectorF dau2(mom2,mom2.massHypothesis(muMass));

    StLorentzVectorF pair;
    pair = dau1+dau2;

    double pt = pair.perp();
    double y = pair.rapidity();
    if(!passMuMuPairCuts(y)) continue;
    double mass = pair.m();
    int dauIsHft1 = muTrk1->isHFTTrack();
    int dauIsHft2 = muTrk2->isHFTTrack();
    int charge1 = muTrk1->charge();
    int charge2 = muTrk2->charge();
    //if(centrality<=0) continue;
    if(charge1!=charge2){ 
      hMuMuNumInvMassvsPtMB->Fill(pt,mass);
      hMuMuNumInvMassvsPt->Fill(pt,mass,centrality);
      if(dauIsHft1&&dauIsHft2){ 
        hMuMuNumInvMassvsPtMBwHft->Fill(pt,mass);
        hMuMuNumInvMassvsPtwHft->Fill(pt,mass,centrality);
      }
    }
    if(charge1==1&&charge2==1){ 
      hMuMuDenInvMassvsPtLikePosMB->Fill(pt,mass);
      hMuMuDenInvMassvsPtLikePos->Fill(pt,mass,centrality);
      if(dauIsHft1&&dauIsHft2){ 
        hMuMuDenInvMassvsPtLikePosMBwHft->Fill(pt,mass);
        hMuMuDenInvMassvsPtLikePoswHft->Fill(pt,mass,centrality);
      }
    }
    if(charge1==-1&&charge2==-1){ 
      hMuMuDenInvMassvsPtLikeNegMB->Fill(pt,mass);
      hMuMuDenInvMassvsPtLikeNeg->Fill(pt,mass,centrality);
      if(dauIsHft1&&dauIsHft2){ 
        hMuMuDenInvMassvsPtLikeNegMBwHft->Fill(pt,mass);
        hMuMuDenInvMassvsPtLikeNegwHft->Fill(pt,mass,centrality);
      }
    }
  }


  hNe->Fill(current_nePlus,current_neMinus);
  hNemu->Fill(current_nePlus+current_neMinus, current_nmuPlus+current_nmuMinus);
  hNmu->Fill(current_nmuPlus, current_nmuMinus);

  makeEEMixedPairs(magBufferPointer,cenBufferPointer, vzBufferPointer,eveBufferPointer);
  makeEMuMixedPairs(magBufferPointer,cenBufferPointer, vzBufferPointer,eveBufferPointer);
  makeMuMuMixedPairs(magBufferPointer,cenBufferPointer, vzBufferPointer,eveBufferPointer);
  copyCurrentToBuffer(magBufferPointer,cenBufferPointer, vzBufferPointer,eveBufferPointer);

  return kStOK;
}//end of main fucntion

bool StMyAnaTreeMaker::passHTEIDCuts(StElectronTrack *eTrk) {
  double pt = eTrk->pMom().perp();
  double p = eTrk->pMom().mag();
  double eta = eTrk->pMom().pseudoRapidity();
  int nHitsFit = eTrk->nHitsFit();
  int nHitsDedx = eTrk->nHitsDedx();
  //int nHitsMax= eTrk->nHitsMax();
  double beta = eTrk->beta();
  double nSigE = eTrk->nSigmaElectron();
  double dca = eTrk->dca();
  double localY = eTrk->localY();
  //double localZ = eTrk->localZ();
  //double ratio = nHitsFit*1./nHitsMax;
  double pve = eTrk->pve();
  double adc0 = eTrk->adc0();
  int nEta = eTrk->nEta();
  int nPhi = eTrk->nPhi();
  double zDist = eTrk->zDist();
  double phiDist = eTrk->phiDist();

  if(pt<mEmcEPtCut[0] || pt>mEmcEPtCut[1]) return false;
  if(eta<mEEtaCut[0] || eta>mEEtaCut[1]) return false;
  //if(dca<mEmcEDcaCut[0] || dca>mEmcEDcaCut[1]) return false;
  if(nHitsFit<mnHitsFitCut[0]||nHitsFit>mnHitsFitCut[1]) return false;
  if(nHitsDedx<mnHitsDedxCut[0]||nHitsDedx>mnHitsDedxCut[1]) return false;
  //if(ratio<mRatioCut[0]||ratio>mRatioCut[1]) return false;

  if(nSigE<mHTEnSigECut[0]||nSigE>mHTEnSigECut[1]) return false;
  if(pve<mEmcEPveCut[0]||pve>mEmcEPveCut[1]) return false;
  //if(nEta<mEnEtaCut[0]||nEta>mEnEtaCut[1]) return false;
  //if(nPhi<mEnPhiCut[0]||nPhi>mEnPhiCut[1]) return false;
  //if(zDist<mEZDistCut[0]||zDist>mEZDistCut[1]) return false;
  //if(phiDist<mEPhiDistCut[0]||phiDist>mEPhiDistCut[1]) return false;

  return true;
}

bool StMyAnaTreeMaker::isHTTrigE(StElectronTrack *eTrk) {
  int emcTrgId = eTrk->emcTriggerId();
  double dsmadc0 = 0;
  double adc0 = eTrk->adc0();
  double pt = eTrk->pMom().perp();
  //if(emcTrgId>=0){
  //	StEmcTrigger *emcTrg = (StEmcTrigger*)mAnaTree->emcTrigger(emcTrgId);
  //	if(!emcTrg) return false;
  //	dsmadc0 = emcTrg->adc0();
  //}else{
  //	return false;
  //}
  if((mTrigSelect==-1||mTrigSelect==1||mTrigSelect==2||mTrigSelect==3||mTrigSelect==4||mTrigSelect==5)&&adc0>mHTAdc0th&&pt>mEmcPtth)  return true;
  else return false;
}


bool StMyAnaTreeMaker::passEIDCuts(StElectronTrack *eTrk) {
  double pt = eTrk->pMom().perp();
  double p = eTrk->pMom().mag();
  double eta = eTrk->pMom().pseudoRapidity();
  int nHitsFit = eTrk->nHitsFit();
  int nHitsDedx = eTrk->nHitsDedx();
  //int nHitsMax= eTrk->nHitsMax();
  double beta = eTrk->beta();
  double nSigE = eTrk->nSigmaElectron();
  double dca = eTrk->dca();
  double localY = eTrk->localY();
  //double localZ = eTrk->localZ();

  if(beta<=0) return false;

  double invBeta = beta!=0 ? 1./beta : 0;
  //double ratio = nHitsFit*1./nHitsMax;

  if(pt<mEPtCut[0] || pt>mEPtCut[1]) return false;
  if(eta<mEEtaCut[0] || eta>mEEtaCut[1]) return false;
  if(dca<mEDcaCut[0] || dca>mEDcaCut[1]) return false;
  if(nHitsFit<mnHitsFitCut[0]||nHitsFit>mnHitsFitCut[1]) return false;
  if(nHitsDedx<mnHitsDedxCut[0]||nHitsDedx>mnHitsDedxCut[1]) return false;
  //if(ratio<mRatioCut[0]||ratio>mRatioCut[1]) return false;
  if(invBeta<mEInvBetaCut[0] || invBeta>mEInvBetaCut[1]) return false;
  if(localY<mELocalYCut[0] || localY>mELocalYCut[1]) return false;
  //if(localZ<mELocalZCut[0] || localZ>mELocalZCut[1]) return false;
  if(p<1.0){
    if(nSigE>mEnSigECut[1]||nSigE<(1.5*(p-0.2)-2.0)) return false;
  }else{
    if(nSigE<mEnSigECut[0]||nSigE>mEnSigECut[1]) return false;
  }

  return true;
}

bool StMyAnaTreeMaker::passMuIDCuts(StMuonTrack *muTrk) {
  double pt = muTrk->pMom().perp();
  double eta = muTrk->pMom().pseudoRapidity();
  int charge = muTrk->charge();
  int nHitsFit = muTrk->nHitsFit();
  int nHitsDedx = muTrk->nHitsDedx();
  double nSigPi = muTrk->nSigmaPion();
  double dca = muTrk->dca();
  double dT = muTrk->deltaTimeOfFlight();
  int channel = muTrk->channel();
  //double dTCor = mtddTCor(dT,channel);
  double dZ = muTrk->deltaZ();
  //double dY = muTrk->deltaY();
  if(pt<mMuPtCut[0] || pt>mMuPtCut[1]) return false;
  if(eta<mMuEtaCut[0] || eta>mMuEtaCut[1]) return false;
  //if(dT<mMudTCut[0] || dT>mMudTCut[1]) return false;
  //if(dTCor<fMudTCutLow->Eval(pt) || dTCor>fMudTCutHigh->Eval(pt)) return false;
  if(dT<fMudTCutLow->Eval(pt) || dT>fMudTCutHigh->Eval(pt)) return false;
  //if(dY<mMudYCut[0] || dY>mMudYCut[1]) return false;
  if(dZ<mMudZCut[0] || dZ>mMudZCut[1]) return false;
  if(dca<mMuDcaCut[0] || dca>mMuDcaCut[1]) return false;
  if(nSigPi<mMunSigPiCut[0] || nSigPi>mMunSigPiCut[1]) return false;
  int triggerFlag = muTrk->triggerFlag();
  if(mTrigSelect==5||mTrigSelect==6||mTrigSelect==7){
    if(triggerFlag>0) return true;
    else return false;
  }	
  return true;
}

bool StMyAnaTreeMaker::passEEPairCuts(double y) {
  if(y<mEEYCut[0] || y>mEEYCut[1]) return false;
  return true;	
}

bool StMyAnaTreeMaker::passEMuPairCuts(double y) {
  if(y<mEMuYCut[0] || y>mEMuYCut[1]) return false;
  return true;	
}

bool StMyAnaTreeMaker::passMuMuPairCuts(double y) {
  if(y<mMuMuYCut[0] || y>mMuMuYCut[1]) return false;
  return true;	
}

int StMyAnaTreeMaker::getCentrality(){
  int gRefMult = mAnaTree->event()->grefMult();	
  int Centrality = 0;
  //temporary
  int cent[] = {10,21,40,71,116,179,263,373,441};
  if(     gRefMult < cent[0]) Centrality = 0;
  else if(gRefMult < cent[1]) Centrality = 1;
  else if(gRefMult < cent[2]) Centrality = 2;
  else if(gRefMult < cent[3]) Centrality = 3;
  else if(gRefMult < cent[4]) Centrality = 4;
  else if(gRefMult < cent[5]) Centrality = 5;
  else if(gRefMult < cent[6]) Centrality = 6;
  else if(gRefMult < cent[7]) Centrality = 7;
  else if(gRefMult < cent[8]) Centrality = 8;
  else Centrality = 9;
  return Centrality;
}
void StMyAnaTreeMaker::printCuts(){

  LOG_INFO<<"analysis cuts:"<<endm;
  LOG_INFO<<mVzCut[0]<<"<mVzCut<"<<mVzCut[1]<<endm;
  LOG_INFO<<mVzDiffCut[0]<<"<mVzDiffCut<"<<mVzDiffCut[1]<<endm;
  LOG_INFO<<"mTrigSelect="<<mTrigSelect<<endm;
  LOG_INFO<<"mHTth="<<mHTth<<endm;

  LOG_INFO<<mnHitsFitCut[0]<<"<mnHitsFitCut<"<<mnHitsFitCut[1]<<endm;
  LOG_INFO<<mnHitsDedxCut[0]<<"<mnHitsDedxCut<"<<mnHitsDedxCut[1]<<endm;
  LOG_INFO<<mRatioCut[0]<<"<mRatioCut<"<<mRatioCut[1]<<endm;
  LOG_INFO<<"mHFTTrackCut="<<mHFTTrackCut<<endm;

  LOG_INFO<<mEPtCut[0]<<"<mEPtCut<"<<mEPtCut[1]<<endm;
  LOG_INFO<<mEEtaCut[0]<<"<mEEtaCut<"<<mEEtaCut[1]<<endm;
  LOG_INFO<<mEDcaCut[0]<<"<mEDcaCut<"<<mEDcaCut[1]<<endm;
  LOG_INFO<<mEInvBetaCut[0]<<"<mEInvBetaCut<"<<mEInvBetaCut[1]<<endm;
  LOG_INFO<<mELocalYCut[0]<<"<mELocalYCut<"<<mELocalYCut[1]<<endm;
  LOG_INFO<<mELocalZCut[0]<<"<mELocalZCut<"<<mELocalZCut[1]<<endm;
  LOG_INFO<<mEnSigECut[0]<<"<mEnSigECut<"<<mEnSigECut[1]<<endm;

  LOG_INFO<<mEmcEPtCut[0]<<"<mEmcEPtCut<"<<mEmcEPtCut[1]<<endm;
  LOG_INFO<<mEmcEEtaCut[0]<<"<mEmcEEtaCut<"<<mEmcEEtaCut[1]<<endm;
  LOG_INFO<<mEmcEPveCut[0]<<"<mEmcEPveCut<"<<mEmcEPveCut[1]<<endm;

  LOG_INFO<<mEnEtaCut[0]<<"<mEnEtaCut<"<<mEnEtaCut[1]<<endm;
  LOG_INFO<<mEnPhiCut[0]<<"<mEnPhiCut<"<<mEnPhiCut[1]<<endm;
  LOG_INFO<<mEZDistCut[0]<<"<mEZDistCut<"<<mEZDistCut[1]<<endm;
  LOG_INFO<<mEPhiDistCut[0]<<"<mEPhiDistCut<"<<mEPhiDistCut[1]<<endm;

  LOG_INFO<<mMuPtCut[0]<<"<mMuPtCut<"<<mMuPtCut[1]<<endm;
  LOG_INFO<<mMuEtaCut[0]<<"<mMuEtaCut<"<<mMuEtaCut[1]<<endm;
  LOG_INFO<<mMuDcaCut[0]<<"<mMuDcaCut<"<<mMuDcaCut[1]<<endm;
  LOG_INFO<<mMunSigPiCut[0]<<"<mMunSigPiCut<"<<mMunSigPiCut[1]<<endm;
  LOG_INFO<<mMudTCut[0]<<"<mMudTCut<"<<mMudTCut[1]<<endm;
  LOG_INFO<<mMudZCut[0]<<"<mMudZCut<"<<mMudZCut[1]<<endm;
  LOG_INFO<<mMudYCut[0]<<"<mMudYCut<"<<mMudYCut[1]<<endm;

  LOG_INFO<<mDauEPtCut[0]<<"<mDauEPtCut<"<<mDauEPtCut[1]<<endm;
  LOG_INFO<<mDauEEtaCut[0]<<"<mDauEEtaCut<"<<mDauEEtaCut[1]<<endm;
  LOG_INFO<<mDauEDcaToVtxCut[0]<<"<mDauEDcaToVtxCut<"<<mDauEDcaToVtxCut[1]<<endm;
  LOG_INFO<<mDauEDcaDistCut[0]<<"<mDauEDcaDistCut<"<<mDauEDcaDistCut[1]<<endm;

  LOG_INFO<<mDauMuPtCut[0]<<"<mDauMuPtCut<"<<mDauMuPtCut[1]<<endm;
  LOG_INFO<<mDauMuEtaCut[0]<<"<mDauMuEtaCut<"<<mDauMuEtaCut[1]<<endm;
  LOG_INFO<<mDauMuDcaToVtxCut[0]<<"<mDauMuDcaToVtxCut<"<<mDauMuDcaToVtxCut[1]<<endm;

  LOG_INFO<<mCosThetaStarCut[0]<<"<mCosThetaStarCut<"<<mCosThetaStarCut[1]<<endm;
  LOG_INFO<<mPointingAngleCut[0]<<"<mPointingAngleCut<"<<mPointingAngleCut[1]<<endm;
  LOG_INFO<<mPairDcaCut[0]<<"<mPairDcaCut<"<<mPairDcaCut[1]<<endm;
  LOG_INFO<<mPairDecayLCut[0]<<"<mPairDecayLCut<"<<mPairDecayLCut[1]<<endm;
  LOG_INFO<<mEEYCut[0]<<"<mEEYCut<"<<mEEYCut[1]<<endm;
  LOG_INFO<<mPairMassCut[0]<<"<mPairMassCut<"<<mPairMassCut[1]<<endm;

  LOG_INFO<<mEEYCut[0]<<"<mEEYCut<"<<mEEYCut[1]<<endm;
  LOG_INFO<<mEMuYCut[0]<<"<mEMuYCut<"<<mEMuYCut[1]<<endm;
  LOG_INFO<<mMuMuYCut[0]<<"<mMuMuYCut<"<<mMuMuYCut[1]<<endm;

}

Double_t StMyAnaTreeMaker::mtddTCor(double dT, int channel){

  int backleg = channel/60+1;
  int module = (channel%60)/12+1;

  return dT + +mMtdT0Corr[backleg-1][module-1];

}

void StMyAnaTreeMaker::makeEEMixedPairs(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer){

  int htFlag = 0;
  if(mTrigSelect==1||mTrigSelect==2||mTrigSelect==3||mTrigSelect==4){
    htFlag = 1;
  }

  TLorentzVector pair(0,0,0,0);
  for(int iBufferEvent=0;iBufferEvent<nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer];iBufferEvent++) {
    //+-------------------------+
    //| current e+  + buffer e- |
    //+-------------------------+
    for(int i=0;i<current_nePlus;i++) {
      for(int j=0;j<buffer_neMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {

        int flag1 = current_ePlusFlag[i];
        int flag2 = buffer_eMinusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_ePlusIsHFT[i];
        int isHft2 = buffer_eMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2/4;
        if(htFlag){
          if(trgFlag1&&trgFlag2) continue; //both fired trigger
          if(!trgFlag1&&!trgFlag2) continue; // no one triggered
          if(flag1==1&&flag2==1) continue; //TPC+TPC

        }
        pair = current_ePlus[i] + buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEEYCut[0]&&pair.Rapidity()<mEEYCut[1]) {
          hEEDenInvMassvsPtMix->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEEDenInvMassvsPtMixwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current e+ loop
    }//end of buffer e- loop

    //+-------------------------+
    //| current e-  + buffer e+ |
    //+-------------------------+
    for(int i=0;i<current_neMinus;i++) {
      for(int j=0;j<buffer_nePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_eMinusFlag[i];
        int flag2 = buffer_ePlusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_eMinusIsHFT[i];
        int isHft2 = buffer_ePlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2/4;
        if(htFlag){
          if(trgFlag1&&trgFlag2) continue; //both fired trigger
          if(!trgFlag1&&!trgFlag2) continue; // no one triggered
          if(flag1==1&&flag2==1) continue; //TPC+TPC
        }

        pair = current_eMinus[i] + buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEEYCut[0]&&pair.Rapidity()<=mEEYCut[1]) {
          hEEDenInvMassvsPtMix->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEEDenInvMassvsPtMixwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current e- loop
    }// end of buffer e+ loop

    //+-------------------------+
    //| current e+  + buffer e+ |
    //+-------------------------+
    // Mixed-likesign
    for(int i=0;i<current_nePlus;i++) {
      for(int j=0;j<buffer_nePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        char flag1 = current_ePlusFlag[i];
        char flag2 = buffer_ePlusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_ePlusIsHFT[i];
        int isHft2 = buffer_ePlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2/4;
        if(htFlag){
          if(trgFlag1&&trgFlag2) continue; //both fired trigger
          if(!trgFlag1&&!trgFlag2) continue; // no one triggered
          if(flag1==1&&flag2==1) continue; //TPC+TPC
        }

        pair = current_ePlus[i] + buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEEYCut[0]&&pair.Rapidity()<=mEEYCut[1]) {
          hEEDenInvMassvsPtMixLikePos->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEEDenInvMassvsPtMixLikePoswHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current e+ loop
    }// end of buffer e+ loop

    //+-------------------------+
    //| current e-  + buffer e- |
    //+-------------------------+
    for(int i=0;i<current_neMinus;i++) {
      for(int j=0;j<buffer_neMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        char flag1 = current_eMinusFlag[i];
        char flag2 = buffer_eMinusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_eMinusIsHFT[i];
        int isHft2 = buffer_eMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2/4;
        if(htFlag){
          if(trgFlag1&&trgFlag2) continue; //both fired trigger
          if(!trgFlag1&&!trgFlag2) continue; // no one triggered
          if(flag1==1&&flag2==1) continue; //TPC+TPC
        }

        pair = current_eMinus[i] + buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEEYCut[0]&&pair.Rapidity()<=mEEYCut[1]) {
          hEEDenInvMassvsPtMixLikeNeg->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEEDenInvMassvsPtMixLikeNegwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current e- loop
    }// end of buffer e- loop

    //+-------------------------+
  }//end of buffer events loop
}

void StMyAnaTreeMaker::makeEMuMixedPairs(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer){
  int emuFlag = 0;
  if(mTrigSelect==5){
    emuFlag = 1;
  }
  int htFlag = 0;
  if(mTrigSelect==1||mTrigSelect==2||mTrigSelect==3||mTrigSelect==4){
    htFlag = 1;
  }

  TLorentzVector pair(0,0,0,0);
  for(int iBufferEvent=0;iBufferEvent<nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer];iBufferEvent++) {
    //+-------------------------+
    //| current e+  + buffer mu-|
    //+-------------------------+
    for(int i=0;i<current_nePlus;i++) {
      for(int j=0;j<buffer_nmuMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_ePlusFlag[i];
        int flag2 = buffer_muMinusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_ePlusIsHFT[i];
        int isHft2 = buffer_muMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2;
        if(emuFlag){
          if(!trgFlag1||!trgFlag2) continue; 
        }
        if(htFlag){
          if(!trgFlag1) continue;
        }

        pair = current_ePlus[i] + buffer_muMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEMuYCut[0]&&pair.Rapidity()<mEMuYCut[1]) {
          hEMuDenInvMassvsPtMix->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEMuDenInvMassvsPtMixwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of buffer mu+ loop
    }//end of current e- loop

    //+-------------------------+
    //| current e-  + buffer mu+ |
    //+-------------------------+
    for(int i=0;i<current_neMinus;i++) {
      for(int j=0;j<buffer_nmuPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_eMinusFlag[i];
        int flag2 = buffer_muPlusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_eMinusIsHFT[i];
        int isHft2 = buffer_muPlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2;
        if(emuFlag){
          if(!trgFlag1||!trgFlag2) continue; 
        }
        if(htFlag){
          if(!trgFlag1) continue;
        }

        pair = current_eMinus[i] + buffer_muPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEMuYCut[0]&&pair.Rapidity()<=mEMuYCut[1]) {
          hEMuDenInvMassvsPtMix->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEMuDenInvMassvsPtMixwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }// end of buffer mu+ loop
    }//end of current e- loop

    //+-------------------------+
    //| current e+  + buffer mu+ |
    //+-------------------------+
    // Mixed-likesign
    for(int i=0;i<current_nePlus;i++) {
      for(int j=0;j<buffer_nmuPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_ePlusFlag[i];
        int flag2 = buffer_muPlusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_ePlusIsHFT[i];
        int isHft2 = buffer_muPlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2;
        if(emuFlag){
          if(!trgFlag1||!trgFlag2) continue; 
        }
        if(htFlag){
          if(!trgFlag1) continue;
        }

        pair = current_ePlus[i] + buffer_muPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEMuYCut[0]&&pair.Rapidity()<=mEMuYCut[1]) {
          hEMuDenInvMassvsPtMixLikePos->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEMuDenInvMassvsPtMixLikePoswHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }// end of buffer mu+ loop
    }//end of current e+ loop

    //+-------------------------+
    //| current e-  + buffer mu- |
    //+-------------------------+
    for(int i=0;i<current_neMinus;i++) {
      for(int j=0;j<buffer_nmuMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_eMinusFlag[i];
        int flag2 = buffer_muMinusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_eMinusIsHFT[i];
        int isHft2 = buffer_muMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1/4;
        int trgFlag2 = flag2;
        if(emuFlag){
          if(!trgFlag1||!trgFlag2) continue; 
        }
        if(htFlag){
          if(!trgFlag1) continue;
        }
        pair = current_eMinus[i] + buffer_muMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mEMuYCut[0]&&pair.Rapidity()<=mEMuYCut[1]) {
          hEMuDenInvMassvsPtMixLikeNeg->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hEMuDenInvMassvsPtMixLikeNegwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }// end of buffer mu- loop
    }//end of current e- loop

  }//end of buffer events loop

}

void StMyAnaTreeMaker::makeMuMuMixedPairs(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer){
  int mumuFlag = 0;
  if(mTrigSelect==5){
    mumuFlag = 1;
  }

  TLorentzVector pair(0,0,0,0);
  for(int iBufferEvent=0;iBufferEvent<nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer];iBufferEvent++) {
    //+-------------------------+
    //| current mu+  + buffer mu-|
    //+-------------------------+
    for(int i=0;i<current_nmuPlus;i++) {
      for(int j=0;j<buffer_nmuMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_muPlusFlag[i];
        int flag2 = buffer_muMinusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_muPlusIsHFT[i];
        int isHft2 = buffer_muMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1;
        int trgFlag2 = flag2;
        if(mumuFlag==1){
          if(!trgFlag1||!trgFlag2) continue;
        }

        pair = current_muPlus[i] + buffer_muMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mMuMuYCut[0]&&pair.Rapidity()<mMuMuYCut[1]) {
          hMuMuDenInvMassvsPtMix->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hMuMuDenInvMassvsPtMixwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current mu+ loop
    }//end of buffer mu- loop

    //+-------------------------+
    //| current mu-  + buffer mu+ |
    //+-------------------------+
    for(int i=0;i<current_nmuMinus;i++) {
      for(int j=0;j<buffer_nmuPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_muMinusFlag[i];
        int flag2 = buffer_muPlusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_muMinusIsHFT[i];
        int isHft2 = buffer_muPlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1;
        int trgFlag2 = flag2;
        if(mumuFlag==1){
          if(!trgFlag1||!trgFlag2) continue;
        }

        pair = current_muMinus[i] + buffer_muPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mMuMuYCut[0]&&pair.Rapidity()<=mMuMuYCut[1]) {
          hMuMuDenInvMassvsPtMix->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hMuMuDenInvMassvsPtMixwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current mu- loop
    }// end of buffer mu+ loop

    //+-------------------------+
    //| current mu+  + buffer mu+ |
    //+-------------------------+
    // Mixed-likesign
    for(int i=0;i<current_nmuPlus;i++) {
      for(int j=0;j<buffer_nmuPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_muPlusFlag[i];
        int flag2 = buffer_muPlusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_muPlusIsHFT[i];
        int isHft2 = buffer_muPlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1;
        int trgFlag2 = flag2;
        if(mumuFlag==1){
          if(!trgFlag1||!trgFlag2) continue;
        }


        pair = current_muPlus[i] + buffer_muPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mMuMuYCut[0]&&pair.Rapidity()<=mMuMuYCut[1]) {
          hMuMuDenInvMassvsPtMixLikePos->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hMuMuDenInvMassvsPtMixLikePoswHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current mu+ loop
    }// end of buffer mu+ loop

    //+-------------------------+
    //| current mu-  + buffer mu- |
    //+-------------------------+
    for(int i=0;i<current_nmuMinus;i++) {
      for(int j=0;j<buffer_nmuMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++) {
        int flag1 = current_muMinusFlag[i];
        int flag2 = buffer_muMinusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int isHft1 = current_muMinusIsHFT[i];
        int isHft2 = buffer_muMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        int trgFlag1 = flag1;
        int trgFlag2 = flag2;
        if(mumuFlag==1){
          if(!trgFlag1||!trgFlag2) continue;
        }


        pair = current_muMinus[i] + buffer_muMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
        if(pair.Rapidity()>mMuMuYCut[0]&&pair.Rapidity()<=mMuMuYCut[1]) {
          hMuMuDenInvMassvsPtMixLikeNeg->Fill(pair.Pt(), pair.M(), current_centrality); 
          if(isHft1&&isHft2) hMuMuDenInvMassvsPtMixLikeNegwHft->Fill(pair.Pt(), pair.M(), current_centrality); 
        }
      }//end of current mu- loop
    }// end of buffer mu- loop

  }//end of buffer events loop
}

void StMyAnaTreeMaker::copyCurrentToBuffer(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer){

  if(nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer]>=nMaxEventsInBuffer) buffer_fullFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer] = kTRUE;

  TRandom3 *gRandom = new TRandom3(iran++);
  int eventPointer = -1;
  if(buffer_fullFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer]){//full, random replace one
    eventPointer =  (int) gRandom->Uniform(0,nMaxEventsInBuffer-1e-6);
  } else {
    eventPointer = nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer];
  }
  delete gRandom;

  buffer_nePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nePlus;
  for(int i=0;i<current_nePlus;i++) {
    buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_ePlus[i];
    buffer_ePlusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_ePlusFlag[i];
    buffer_ePlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_ePlusIsHFT[i];
  }

  buffer_neMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_neMinus;
  for(int i=0;i<current_neMinus;i++) {
    buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_eMinus[i];
    buffer_eMinusFlag[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_eMinusFlag[i];
    buffer_eMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_eMinusIsHFT[i];
  }

  buffer_nmuPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nmuPlus;
  for(int i=0;i<current_nmuPlus;i++) {
    buffer_muPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_muPlus[i];
    buffer_muPlusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_muPlusIsHFT[i];
  }

  buffer_nmuMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nmuMinus;
  for(int i=0;i<current_nmuMinus;i++) {
    buffer_muMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_muMinus[i];
    buffer_muMinusIsHFT[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_muMinusIsHFT[i];
  }

  if(nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer]<nMaxEventsInBuffer) {
    nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer][eveBufferPointer] += 1;
  }
}

bool StMyAnaTreeMaker::passHadronCuts(StHadronTrack* hTrk)
{
  int charge = hTrk->charge();
  double pt = hTrk->pMom().perp();
  double eta = hTrk->pMom().pseudoRapidity();
  double phi = hTrk->pMom().phi();
  double dca = hTrk->dca();

  if(charge!=0 && 
      mHadPtCut[0] < pt  && pt  < mHadPtCut[1] &&
      mHadEtaCut[0]< eta && eta < mHadEtaCut[1] &&
      mHadDcaCut[0]< dca && dca < mHadDcaCut[1]) 
  {
    return true;
  }
}

double StMyAnaTreeMaker::delPhiCorrect(double p)
{
  double PI = 3.1415926;  
  if(p < -PI/2.0) p+=2*PI;
  if(p >  3.*PI/2.)  p-=2.*PI;
  return p;
}

