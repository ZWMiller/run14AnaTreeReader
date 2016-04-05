#include "StCutsTreeMaker.h"
#include "StRoot/StPicoAnaTreeMaker/StAnaTree.h"
#include "StRoot/StPicoAnaTreeMaker/StPicoAnaTreeMaker.h"
#include "StRoot/StPicoAnaTreeMaker/StEventHeader.h"
#include "StRoot/StPicoAnaTreeMaker/StElectronTrack.h"
//#include "StRoot/StPicoAnaTreeMaker/StMuonTrack.h"
#include "StRoot/StPicoAnaTreeMaker/StEEPair.h"
//#include "StRoot/StPicoAnaTreeMaker/StEMuPair.h"
//#include "StRoot/StPicoAnaTreeMaker/StMuMuPair.h"
#include "StRoot/StPicoAnaTreeMaker/StEmcTrigger.h"

#include "StMuDSTMaker/COMMON/StMuEmcUtil.h"
#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StEmcUtil/geometry/StEmcGeom.h"

#include "StLorentzVector.hh"
#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "mBadRunList.h"
#include "algorithm"
#include "TRandom3.h"


#define nCentrals 10
#define eMass 0.000510999

#define nMassBins 500
#define massMin 0
#define massMax 5 
#define nPtBins 200
#define ptMin 0
#define ptMax 12 


ClassImp(StCutsTreeMaker)

	//-----------------------------------------------------------------------------
	StCutsTreeMaker::StCutsTreeMaker(const char* name, StPicoAnaTreeMaker *treeMaker, const char* outName)
: StMaker(name)
{
	mPicoAnaTreeMaker = treeMaker;
	mAnaTree = 0;
	TH1F:: SetDefaultSumw2();//zaochen add
	mOutName = outName;

	mTrigSelect = 0;//trigger selection code, 
	if(mTrigSelect==0){
		mVzCut[0] = -6; mVzCut[1] = 6;
		mVzDiffCut[0] = -3; mVzDiffCut[1] = 3;
	}
	else if(mTrigSelect==2||mTrigSelect==3){
		mVzCut[0] = -30; 
		mVzCut[1] = 30;
		mVzDiffCut[0] = -4; mVzDiffCut[1] = 4;
	}else{
		mVzCut[0] = -200; 
		mVzCut[1] = 200;
		mVzDiffCut[0] = -4; mVzDiffCut[1] = 4;
	}
	mVxCut[0] = -1.0e-5; mVxCut[1] =1.0e-5;
	mVyCut[0] = -1.0e-5; mVyCut[1] =1.0e-5;

	mnHitsFitCut[0] = 25; mnHitsFitCut[1] = 50;
	mnHitsDedxCut[0] = 15; mnHitsDedxCut[1] = 50;
	mRatioCut[0] = 0.52; mRatioCut[1] = 1.;

	mEPtCut[0] = 0.2;             mEPtCut[1] = 30;
	mEPCut[0] = 0.2;              mEPCut[1] = 30;

	mEEtaCut[0] = -1.;            mEEtaCut[1] = 1.;
	mEDcaCut[0] = 0.;             mEDcaCut[1] = 3.;
	mEInvBetaCut[0] = 0.97;       mEInvBetaCut[1] = 1.03;
	mELocalYCut[0] = -1.8;        mELocalYCut[1] = 1.8;
	mELocalZCut[0] = -3.05;       mELocalZCut[1] = 3.05;
	mEnSigECut[0] = -1.5;         mEnSigECut[1] = 2.5;

	mPairMassCut[0]=2.0;          mPairMassCut[1]=5.0;  
	mPEMassCut[0] = 0.;           mPEMassCut[1] = 0.1;
	mEEPairYCut[0]=-1.0;           mEEPairYCut[1]=1.0;
	mNBadRuns = sizeof(mBadRuns)/sizeof(int);



}

//----------------------------------------------------------------------------- 
StCutsTreeMaker::~StCutsTreeMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StCutsTreeMaker::Init() {

	if(mOutName!="") {
		fout = new TFile(mOutName.Data(),"RECREATE");
	}else{
		fout = new TFile("test.ana.root","RECREATE");
	}
	declareTree();
	declareHistograms();


	return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StCutsTreeMaker::Finish() {
	fout->cd();
	fout->Write();
	fout->Close();
	printCuts();
	return kStOK;
}


//-----------------------------------------------------------------------------
void StCutsTreeMaker::declareTree() {

	fout->cd();


	mb_Tree = new TTree("mb_Tree","mb_Tree");
	mb_Tree->SetAutoSave(100000);

	//Event level information                                                             
	//define branches                                                                     
	mb_Tree->Branch( "mb_RefMult",  &mb_RefMult,  "mb_RefMult/I" );
	mb_Tree->Branch( "mb_GRefMult", &mb_GRefMult, "mb_GRefMult/I" );
	mb_Tree->Branch( "mb_VzTPC",    &mb_VzTPC,    "mb_VzTPC/F" );
	mb_Tree->Branch( "mb_Centrality", &mb_Centrality,  "mb_Centrality/F" );
	mb_Tree->Branch( "mb_Ranking",    &mb_Ranking,     "mb_Ranking/F" );

	//electron pairs
	mb_Tree->Branch( "mb_EEPairSign",  &mb_EEPairSign,  "mb_EEPairSign/I" );
	mb_Tree->Branch( "mb_EEPairPETag",  &mb_EEPairPETag,  "mb_EEPairPETag/I" );

	mb_Tree->Branch( "mb_EEPairgPt",  &mb_EEPairgPt,  "mb_EEPairgPt/F" );
	mb_Tree->Branch( "mb_EEPairgPz",  &mb_EEPairgPz,  "mb_EEPairgPz/F" );
	mb_Tree->Branch( "mb_EEPairgEta",  &mb_EEPairgEta,  "mb_EEPairgEta/F" );
	mb_Tree->Branch( "mb_EEPairgPhi",  &mb_EEPairgPhi,  "mb_EEPairgPhi/F" );
	mb_Tree->Branch( "mb_EEPairgMass", &mb_EEPairgMass,    "mb_EEPairgMass/F" );
	mb_Tree->Branch( "mb_EEPairgY",   &mb_EEPairgY,   "mb_EEPairgY/F" );

	mb_Tree->Branch( "mb_EEPairPPt",  &mb_EEPairPPt,  "mb_EEPairPPt/F" );
	mb_Tree->Branch( "mb_EEPairPPz",  &mb_EEPairPPz,  "mb_EEPairPPz/F" );
	mb_Tree->Branch( "mb_EEPairPEta",  &mb_EEPairPEta,  "mb_EEPairPEta/F" );
	mb_Tree->Branch( "mb_EEPairPPhi",  &mb_EEPairPPhi,  "mb_EEPairPPhi/F" );
	mb_Tree->Branch( "mb_EEPairPMass", &mb_EEPairPMass,    "mb_EEPairPMass/F" );
	mb_Tree->Branch( "mb_EEPairPY",   &mb_EEPairPY,   "mb_EEPairPY/F" );


	mb_Tree->Branch( "mb_EEPairOx",   &mb_EEPairOx,   "mb_EEPairOx/F" );
	mb_Tree->Branch( "mb_EEPairOy",   &mb_EEPairOy,   "mb_EEPairOy/F" );
	mb_Tree->Branch( "mb_EEPairOz",   &mb_EEPairOz,   "mb_EEPairOz/F" );
	mb_Tree->Branch( "mb_EEPairOr",   &mb_EEPairOr,   "mb_EEPairOr/F" );


	mb_Tree->Branch( "mb_EEPairEEDca",      &mb_EEPairEEDca,    "mb_EEPairEEDca/F" );
	mb_Tree->Branch( "mb_EEPairDecayLth",   &mb_EEPairDecayLth, "mb_EEPairDecayLth/F" );
	mb_Tree->Branch( "mb_EEPairDca",        &mb_EEPairDca,      "mb_EEPairDca/F" );
	mb_Tree->Branch( "mb_EEPairCtau",       &mb_EEPairCtau,     "mb_EEPairCtau/F" );
	mb_Tree->Branch( "mb_EEPairCosThetaStar",  &mb_EEPairCosThetaStar, "mb_EEPairCosThetaStar/F" );
	mb_Tree->Branch( "mb_EEPairPTAngle",       &mb_EEPairPTAngle,      "mb_EEPairPTAngle/F" );
	mb_Tree->Branch( "mb_EEPairDPAngle",       &mb_EEPairDPAngle,      "mb_EEPairDPAngle/F" );
	mb_Tree->Branch( "mb_EEPairOpenAngle",       &mb_EEPairOpenAngle,      "mb_EEPairOpenAngle/F" );
	mb_Tree->Branch( "mb_EEPairPhiV",       &mb_EEPairPhiV,      "mb_EEPairPhiV/F" );

	//3. daughter electron 
	mb_Tree->Branch( "mb_DauEisHFT",    &mb_DauEisHFT,   "mb_DauEisHFT[2]/O");
	mb_Tree->Branch( "mb_DauEisTpc",    &mb_DauEisTpc,   "mb_DauEisTpc[2]/O");
	mb_Tree->Branch( "mb_DauEisEmc",    &mb_DauEisEmc,   "mb_DauEisEmc[2]/O");
	mb_Tree->Branch( "mb_DauEisTrig",    &mb_DauEisTrig,   "mb_DauEisTrig[2]/O");

	mb_Tree->Branch( "mb_DauEgPt", &mb_DauEgPt,  "mb_DauEgPt[2]/F" );
	mb_Tree->Branch( "mb_DauEgPz", &mb_DauEgPz,  "mb_DauEgPz[2]/F" );
	mb_Tree->Branch( "mb_DauEgEta", &mb_DauEgEta,  "mb_DauEgEta[2]/F" );
	mb_Tree->Branch( "mb_DauEgPhi", &mb_DauEgPhi,  "mb_DauEgPhi[2]/F" );

	mb_Tree->Branch( "mb_DauEPt",   &mb_DauEPt,  "mb_DauEPt[2]/F" );
	mb_Tree->Branch( "mb_DauEPz",   &mb_DauEPz,  "mb_DauEPz[2]/F" );
	mb_Tree->Branch( "mb_DauEpEta",  &mb_DauEpEta,  "mb_DauEpEta[2]/F" );
	mb_Tree->Branch( "mb_DauEpPhi",  &mb_DauEpPhi,  "mb_DauEpPhi[2]/F" );

	mb_Tree->Branch( "mb_DauEQ",         &mb_DauEQ,        "mb_DauEQ[2]/I");
	mb_Tree->Branch( "mb_DauEDcaXY",     &mb_DauEDcaXY,    "mb_DauEDcaXY[2]/F");
	mb_Tree->Branch( "mb_DauEDcaZ",      &mb_DauEDcaZ,     "mb_DauEDcaZ[2]/F");

	mb_Tree->Branch( "mb_DauEnSigE",     &mb_DauEnSigE,    "mb_DauEnSigE[2]/F");
	mb_Tree->Branch( "mb_DauEDca",       &mb_DauEDca,      "mb_DauEDca[2]/F");
	mb_Tree->Branch( "mb_DauEInvBeta",   &mb_DauEInvBeta,    "mb_DauEInvBeta[2]/F");
	mb_Tree->Branch( "mb_DauELocalY",    &mb_DauELocalY,     "mb_DauELocalY[2]/F");

	mb_Tree->Branch( "mb_DauEAdc0",     &mb_DauEAdc0,    "mb_DauEAdc0[2]/I" );
	mb_Tree->Branch( "mb_DauEe0",       &mb_DauEe0,      "mb_DauEe0[2]/F");
	mb_Tree->Branch( "mb_DauEe",        &mb_DauEe,       "mb_DauEe[2]/F");
	mb_Tree->Branch( "mb_DauEpve",      &mb_DauEpve,     "mb_DauEpve[2]/F");
	mb_Tree->Branch( "mb_DauEnEta",     &mb_DauEnEta,    "mb_DauEnEta[2]/I");
	mb_Tree->Branch( "mb_DauEnPhi",     &mb_DauEnPhi,    "mb_DauEnPhi[2]/I");
	mb_Tree->Branch( "mb_DauEzDist",    &mb_DauEzDist,   "mb_DauEzDist[2]/F");
	mb_Tree->Branch( "mb_DauEphiDist",  &mb_DauEphiDist, "mb_DauEphiDist[2]/F");
	//----------------------------------------------------------------------
}
//--------------------------

void StCutsTreeMaker::declareHistograms() {

	hnEvents = new TH1F("hnEvents","hnEvents",10,0,10);
	hnTracks = new TH1F("hnTracks","hnTracks",10,0,10);   
	hNe = new TH2F("hNe","#e+ vs. #e-;#e^{+} candidate;#e^{-} candidate;Counts",100,0,100,100,0,100);
	hEventPlane0 = new TH1F("hEventPlane0", "original eventplane;original #Psi;", 200, 0.0, 3.15);
	hEventPlane_ReCen = new TH1F("hEventPlane_ReCen", "After recenter correction; ReCen #Psi;", 200, 0.0, 3.15);


}

//-----------------------------------------------------------------------------
//----------------------------------------------------------------------------- 
void StCutsTreeMaker::Clear(Option_t *opt) {

}

//----------------------------------------------------------------------------- 
Int_t StCutsTreeMaker::Make() {

	if(!mPicoAnaTreeMaker) {
		LOG_WARN << " No PicoAnaTreeMaker! Skip! " << endm;
		return kStOK;
	}

	mAnaTree = mPicoAnaTreeMaker->anaTree();

	if(!mAnaTree) {
		LOG_WARN << " No AnaTree! Skip! " << endm;
		return kStOK;
	}


	if(!mAnaTree->event()){
		//LOG_WARN << " No event! Skip! " << endm;
		return kStOK;
	}
	hnEvents->Fill(0);
	//	StEmcGeom *mEmcGeom = StEmcGeom::instance("bemc");

	//1. reject bad run
	int runId = mAnaTree->event()->runId();
	for(int i=0;i<mNBadRuns;i++){
		if(runId==mBadRuns[i]) return kStOK;
	}

	hnEvents->Fill(1);


	Double_t vzVpd=mAnaTree->event()->vzVpd();
	StThreeVectorF mPrimaryVertex = mAnaTree->event()->primaryVertex();

	//3. vertex cut
	double vx,vy,vz;
	vx=mPrimaryVertex.x();
	vy=mPrimaryVertex.y();
	vz=mPrimaryVertex.z();

	if(mPrimaryVertex.z()<=mVzCut[0]||mPrimaryVertex.z()>=mVzCut[1]) return kStOK; //TPC VZ cut
	if(fabs(vx)<=mVxCut[1]&&fabs(vy)<=mVyCut[1]&&fabs(vz)<=mVxCut[1])return kStOK;// not at zero position
	hnEvents->Fill(2);
	if(vzVpd-mPrimaryVertex.z()<=mVzDiffCut[0]||vzVpd-mPrimaryVertex.z()>=mVzDiffCut[1]) return kStOK;//dVzCut

	hnEvents->Fill(3);
	//1. fill event level tree information
	mb_RefMult=mAnaTree->event()->refMult();
	mb_GRefMult=mAnaTree->event()->grefMult();

	mb_Ranking=mAnaTree->event()->ranking();
	mb_VzTPC=vz;

	int centrality = getCentrality();
	mb_Centrality=centrality;





	//add for the eventplane shift factor

	float evtPlane0 = mAnaTree->event()->eventplane0();
	float evtPlane_ReCen = mAnaTree->event()->eventPlane();


	hEventPlane0->Fill(evtPlane0);
	hEventPlane_ReCen->Fill(evtPlane_ReCen);

	for(int k=1; k<21; k++){

		shiftcorrectioncos[k-1] += cos( 2*k*evtPlane_ReCen );
		shiftcorrectionsin[k-1] += sin( 2*k*evtPlane_ReCen );


	}

//	for(int i=0; i<20; i++){
//		hProShiftCorr_Cos[i]->Fill( shiftcorrectioncos[i] );
//		hProShiftCorr_Sin[i]->Fill( shiftcorrectionsin[i] );
//
//
//	}
	//------------------------------------------















/*

	//2.loop over EEpairs 
	int nEEPairs = mAnaTree->numberOfEEPairs();
	//----------------------------------------
	//loop to tag the photonic electron pairs                      
	vector<int> PEindexLIST;
	for(int iee=0; iee< nEEPairs;iee++){
		StEEPair* eepair = (StEEPair*) mAnaTree->eePair(iee);
		int dauIndex1 = eepair->dauIndex1();
		int dauIndex2 = eepair ->dauIndex2();

		int charge1 = eepair->dauCharge1();
		int charge2 = eepair->dauCharge2();

		double mass = eepair->pairPMass();

		if(mass>mPEMassCut[0]&&mass<mPEMassCut[1]){
			if(charge1!=charge2){
				PEindexLIST.push_back(dauIndex1);
				PEindexLIST.push_back(dauIndex2);
				// cout<<"dauIndex1/dauIndex2:"<<dauIndex1<<"/"<<dauIndex2<<endl;          
			}
		}
	}



	for(int iee=0; iee< nEEPairs;iee++){	
		StEEPair* ee = (StEEPair*) mAnaTree->eePair(iee);
		int dauIndex1 = ee->dauIndex1();
		int dauIndex2 = ee->dauIndex2();
		StElectronTrack *eTrk1 = mAnaTree->eTrack(dauIndex1);
		StElectronTrack *eTrk2 = mAnaTree->eTrack(dauIndex2);

		//------------------------
		std::vector<int>::iterator position1 = std::find(PEindexLIST.begin(), PEindexLIST.end(), dauIndex1);
		bool exists1 = ( position1 != PEindexLIST.end() );

		std::vector<int>::iterator position2 = std::find(PEindexLIST.begin(), PEindexLIST.end(), dauIndex2);
		bool exists2 = ( position2 != PEindexLIST.end() );
		mb_EEPairPETag=-99;
		if(!exists1&&!exists2) mb_EEPairPETag=0;// only 0 photonic e in this pair
		if((!exists1&&exists2) || (exists1&&!exists2) )mb_EEPairPETag=1; // 1 photonic e in this pair
		if( exists1&&exists2 ) mb_EEPairPETag=2; // 2 photonic e in this pair


		mb_DauEisPE[0]=exists1;
		mb_DauEisPE[1]=exists2;
		//---------------------------------
		int flag = 0;
		if(mTrigSelect==0){
			if(passTPCEIDCuts(eTrk1)&&passTPCEIDCuts(eTrk2)) flag=1;
		}
		if(flag==0) continue;
		if(!passEEPairCuts(ee)) continue;		

		double y = ee->pairPY();
		double pt = ee->pairPPt();
		double mass = ee->pairPMass();
		if(mass<mPairMassCut[0]||mass>mPairMassCut[1]) continue; // select mass range
		//fill pair information


		int charge1 = ee->dauCharge1();
		int charge2 = ee->dauCharge2();
		if(fabs(charge1)!=1)continue;
		if(fabs(charge2)!=1)continue;
		mb_EEPairSign=charge1+charge2;


		StThreeVectorF gmom = ee->pairMom();
		Double_t gpz=gmom.z();
		mb_EEPairgPt= ee->pairPt();
		mb_EEPairgPz= gpz;
		mb_EEPairgY= ee->pairY();
		mb_EEPairgEta= ee->pairEta();
		mb_EEPairgPhi= ee->pairPhi();
		mb_EEPairgMass= ee->pairMass();

		//fill pair information
		StThreeVectorF pmom = ee->pairPMom();
		Double_t ppz=pmom.z();
		mb_EEPairPPt=pt;
		mb_EEPairPPz= ppz;
		mb_EEPairPY= y;
		mb_EEPairPEta= ee->pairPEta();
		mb_EEPairPPhi= ee->pairPPhi();
		mb_EEPairPMass= mass;

		StThreeVectorF origin = ee->pairOrigin() - mPrimaryVertex;
		mb_EEPairOx  =  origin.x();
		mb_EEPairOy  =  origin.y();
		mb_EEPairOz  =  origin.z();
		mb_EEPairOr  =  sqrt(mb_EEPairOx*mb_EEPairOx+mb_EEPairOy*mb_EEPairOy);

		mb_EEPairEEDca =  ee->dauDcaDist();
		mb_EEPairDecayLth = ee->pairDecayL();
		mb_EEPairDca = ee->pairDca();
		mb_EEPairCtau = ee->pairCtau();
		mb_EEPairCosThetaStar =ee->cosThetaStar();
		mb_EEPairPTAngle = ee->pointingAngle();
		mb_EEPairPhiV = ee->pairPhiV();

		double gpt1 = eTrk1->gMom().perp();
		double gpz1 = eTrk1->gMom().z();
		double gpx1 = eTrk1->gMom().x();
		double gpy1 = eTrk1->gMom().y();
		double gpt2 = eTrk2->gMom().perp();
		double gpz2 = eTrk2->gMom().z();
		double gpx2 = eTrk2->gMom().x();
		double gpy2 = eTrk2->gMom().y();
		mb_EEPairDPAngle=-99;
		Double_t tempcos=(gpz1*gpz2+gpt1*gpt2)/sqrt( (gpt1*gpt1+gpz1*gpz1)*(gpt2*gpt2+gpz2*gpz2 ) );
		Double_t tempopencos=( gpx1*gpx2+gpy1*gpy2 )/(gpt1*gpt2);
		mb_EEPairDPAngle = TMath::ACos(tempcos);
		mb_EEPairOpenAngle = TMath::ACos(tempopencos);


		mb_DauEAdc0[0]= eTrk1->adc0();
		mb_DauEAdc0[1]= eTrk2->adc0();

		mb_DauEgPt[0]=gpt1; mb_DauEgPt[1]=gpt2;
		mb_DauEgPz[0]=gpz1; mb_DauEgPz[1]=gpz2;

		mb_DauEPt[0]=eTrk1->pMom().perp();
		mb_DauEPt[1]=eTrk2->pMom().perp();
		mb_DauEPz[0]=eTrk1->pMom().z();
		mb_DauEPz[1]=eTrk2->pMom().z();

		mb_DauEgEta[0]=eTrk1->gEta();   mb_DauEgEta[1]=eTrk2->gEta();
		mb_DauEgPhi[0]=eTrk1->gPhi();   mb_DauEgPhi[1]=eTrk2->gPhi();

		mb_DauEpEta[0]=eTrk1->pMom().pseudoRapidity();   mb_DauEpEta[1]=eTrk2->pMom().pseudoRapidity();
		mb_DauEpPhi[0]=eTrk1->pMom().phi();   mb_DauEpPhi[1]=eTrk2->pMom().phi();


		mb_DauEQ[0]=charge1; mb_DauEQ[1]=charge2;
		mb_DauEnSigE[0]=eTrk1->nSigmaElectron();
		mb_DauEnSigE[1]=eTrk2->nSigmaElectron();
		mb_DauEDca[0]=eTrk1->dca();     mb_DauEDca[1]=eTrk2->dca();
		mb_DauEDcaXY[0]=eTrk1->dcaXY(); mb_DauEDcaXY[1]=eTrk2->dcaXY();
		mb_DauEDcaZ[0]=eTrk1->dcaZ();   mb_DauEDcaZ[1]=eTrk2->dcaZ();

		mb_DauEpve[0]=eTrk1->pve();   mb_DauEpve[1]=eTrk2->pve();
		mb_DauEnEta[0]=eTrk1->nEta(); mb_DauEnEta[1]=eTrk2->nEta();
		mb_DauEnPhi[0]=eTrk1->nPhi(); mb_DauEnPhi[1]=eTrk2->nPhi();

		mb_DauEphiDist[0]=eTrk1->phiDist(); 
		mb_DauEphiDist[1]=eTrk2->phiDist();
		mb_DauEzDist[0]=eTrk1->zDist(); mb_DauEzDist[1]=eTrk2->zDist();

		mb_DauEisHFT[0] = eTrk1->isHFTTrack();
		mb_DauEisHFT[1] = eTrk2->isHFTTrack();


		mb_DauEisTpc[0]=passTPCEIDCuts(eTrk1);
		mb_DauEisTpc[1]=passTPCEIDCuts(eTrk2);

		//		mb_DauEisEmc[0]=passEMCEIDCuts(eTrk1);
		//		mb_DauEisEmc[1]=passEMCEIDCuts(eTrk2);
		//
		//		mb_DauEisTrig[0]=isHTTrigE(eTrk1);
		//		mb_DauEisTrig[1]=isHTTrigE(eTrk2);

		mb_DauEe[0]=eTrk1->e();mb_DauEe[1]=eTrk2->e();
		mb_DauEe0[0]=eTrk1->e0();mb_DauEe0[1]=eTrk2->e0();

		double beta1= eTrk1->beta();
		double invBeta1 = beta1!=0 ? 1./beta1 : 0;	
		double beta2= eTrk2->beta();
		double invBeta2 = beta2!=0 ? 1./beta2 : 0;

		mb_DauEInvBeta[0]=invBeta1; mb_DauEInvBeta[1]=invBeta2;
		mb_DauELocalY[0]=eTrk1->localY();
		mb_DauELocalY[1]=eTrk2->localY();

		mb_Tree->Fill();

	}

*/




	return kStOK;
}//end of main fucntion
//=================================================================

//-------------------------------------------------------------
//===========================================================


//-----------------------------------------                                              

bool StCutsTreeMaker::passTPCEIDCuts(StElectronTrack *eTrk) {
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
	//	double localZ = eTrk->localZ();

	if(beta<=0) return false;

	double invBeta = beta!=0 ? 1./beta : 0;
	//double ratio = (nHitsFit*1.)/(1.0*nHitsMax);


	if(nHitsFit<mnHitsFitCut[0]||nHitsFit>mnHitsFitCut[1]) return false;
	if(nHitsDedx<mnHitsDedxCut[0]||nHitsDedx>mnHitsDedxCut[1]) return false;
	//if(ratio<mRatioCut[0]||ratio>mRatioCut[1]) return false;

	if(pt<mEPtCut[0] || pt>mEPtCut[1]) return false;
	if(eta<mEEtaCut[0] || eta>mEEtaCut[1]) return false;
	if(dca<mEDcaCut[0] || dca>mEDcaCut[1]) return false;
	if(invBeta<mEInvBetaCut[0] || invBeta>mEInvBetaCut[1]) return false;
	if(localY<mELocalYCut[0] || localY>mELocalYCut[1]) return false;

	//if(localZ<mELocalZCut[0] || localZ>mELocalZCut[1]) return false;
	//if( nSigE>mEnSigECut[1]||nSigE<mEnSigECut[0] ) return false;
	if(p<1.0){
		if(nSigE>mEnSigECut[1]||nSigE<(1.5*(p-0.2)-2.0)) return false;
	}else{
		if(nSigE<mEnSigECut[0]||nSigE>mEnSigECut[1]) return false;
	}

	return true;
}
//-------------------------------------------------------------
//-------------------------------------------------------------
bool StCutsTreeMaker::passEEPairCuts(StEEPair *ee) {
	double y = ee->pairY();
	if(y<mEEPairYCut[0] || y>mEEPairYCut[1]) return false;
	return true;	
}
//-------------------------------------------------------
int StCutsTreeMaker::getCentrality(){
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
//-------------------------------------------------------

//==============================================================================

void StCutsTreeMaker::printCuts(){

	LOG_INFO<<"Cuts for mini Tree:"<<endm;
	LOG_INFO<<mVzCut[0]<<"<mVzCut<"<<mVzCut[1]<<endm;
	LOG_INFO<<mVxCut[0]<<"<mVxCut<"<<mVxCut[1]<<endm;
	LOG_INFO<<mVyCut[0]<<"<mVyCut<"<<mVyCut[1]<<endm;
	LOG_INFO<<mVzDiffCut[0]<<"<mVzDiffCut<"<<mVzDiffCut[1]<<endm;

	LOG_INFO<<"mTrigSelect="<<mTrigSelect<<endm;
	LOG_INFO<<mnHitsFitCut[0]<<"<mnHitsFitCut<"<<mnHitsFitCut[1]<<endm;
	LOG_INFO<<mnHitsDedxCut[0]<<"<mnHitsDedxCut<"<<mnHitsDedxCut[1]<<endm;
	LOG_INFO<<mRatioCut[0]<<"<mRatioCut<"<<mRatioCut[1]<<endm;

	LOG_INFO<<mEPtCut[0]<<"<mEPtCut<"<<mEPtCut[1]<<endm;
	LOG_INFO<<mEEtaCut[0]<<"<mEEtaCut<"<<mEEtaCut[1]<<endm;
	LOG_INFO<<mEDcaCut[0]<<"<mEDcaCut<"<<mEDcaCut[1]<<endm;
	LOG_INFO<<mEInvBetaCut[0]<<"<mEInvBetaCut<"<<mEInvBetaCut[1]<<endm;
	LOG_INFO<<mELocalYCut[0]<<"<mELocalYCut<"<<mELocalYCut[1]<<endm;
	LOG_INFO<<mEnSigECut[0]<<"<mEnSigECut<"<<mEnSigECut[1]<<endm;
	LOG_INFO<<mPairMassCut[0]<<"<mPairMassCut<"<<mPairMassCut[1]<<endm;
	LOG_INFO<<mPEMassCut[0]<<"<mEmcEnSigECut<"<<mPEMassCut[1]<<endm;
	LOG_INFO<<mEEPairYCut[0]<<"<mEEPairYCut<"<<mEEPairYCut[1]<<endm;



}
//------------------------------------------------------------------------------
//

