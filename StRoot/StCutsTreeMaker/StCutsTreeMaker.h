#ifndef StCutsTreeMaker_h
#define StCutsTreeMaker_h
//#include "StRoot/StPicoAnaTreeMaker/StPicoAnaTreeMaker.h"
#include "StMaker.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "StPhysicalHelixD.hh"

class StAnaTree;
class StEventHeader;
class StElectronTrack;
//class StMuonTrack;
class StEEPair;
//class StEMuPair;
//class StMuMuPair;
class StPicoAnaTreeMaker;

class TString;
class TH1F;
class TH2F;
class TH3F;
class TFile;
class TF1;
class StCutsTreeMaker : public StMaker {
	public:
		StCutsTreeMaker(const char *name, StPicoAnaTreeMaker *treeMaker, const char *outName);
		virtual ~StCutsTreeMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

		void    declareTree();
		void    declareHistograms();
		bool	passTPCEIDCuts(StElectronTrack *);
		bool	passEEPairCuts(StEEPair *);

		int		getCentrality();
		void	printCuts();
	private:

		StPicoAnaTreeMaker *mPicoAnaTreeMaker;
		StAnaTree          *mAnaTree;

		TString    mOutName;

		TFile*	   fout;

		Int_t mNBadRuns;
		Int_t       mTrigSelect; //-1 - all, 0 - MB, 1 - HT0, 2 - HT1, 3 - HT2, 4 - HT3, 5 - EMu, 6 - dimuon.. 		
		TH1F *hnEvents;
		TH1F *hnTracks;
		TH2F *hNe;



		TH1F *hEventPlane0;
		TH1F *hEventPlane_ReCen;
		float shiftcorrectioncos[20], shiftcorrectionsin[20];



		//----------------------------------------------------------------------------------
		TTree*     mb_Tree;//store the same event pairs, the real pairs
		//Event level information 
		Int_t         mb_RefMult; 
		Int_t         mb_GRefMult; 
		Float_t       mb_VzTPC;
		Float_t       mb_Ranking; 
		Float_t       mb_Centrality;

		//ee pairs information
		Int_t         mb_EEPairPETag;
		Int_t         mb_EEPairSign;  //-1=unlike sign , +1=likesign
		//1.1 reconstructed from global tracks
		Float_t       mb_EEPairgPt;
		Float_t       mb_EEPairgPz;
		Float_t       mb_EEPairgY;
		Float_t       mb_EEPairgEta;
		Float_t       mb_EEPairgPhi;
		Float_t       mb_EEPairgMass;

		//1.2 reconstructed from primary tracks
		Float_t       mb_EEPairPPt;
		Float_t       mb_EEPairPPz;
		Float_t       mb_EEPairPEta;
		Float_t       mb_EEPairPPhi;
		Float_t       mb_EEPairPMass;
		Float_t       mb_EEPairPY;

		//2.
		Float_t       mb_EEPairOx;
		Float_t       mb_EEPairOy;
		Float_t       mb_EEPairOz;
		Float_t       mb_EEPairOr;

		Float_t       mb_EEPairEEDca;
		Float_t       mb_EEPairDecayLth;
		Float_t       mb_EEPairDca;
		Float_t       mb_EEPairCtau;
		Float_t       mb_EEPairCosThetaStar;
		Float_t       mb_EEPairPTAngle;
		Float_t       mb_EEPairDPAngle;
		Float_t       mb_EEPairOpenAngle;
		Float_t       mb_EEPairPhiV;

		//another angle Bingchu mentioned
		//3. daughter electron

		Bool_t        mb_DauEisTpc[2];
		Bool_t        mb_DauEisEmc[2];
		Bool_t        mb_DauEisTrig[2];
		Bool_t        mb_DauEisPE[2];
		Bool_t        mb_DauEisHFT[2];

		Int_t         mb_DauEAdc0[2];

		Float_t       mb_DauEgPt[2];
		Float_t       mb_DauEgPz[2];
		Float_t       mb_DauEgEta[2];
		Float_t       mb_DauEgPhi[2];

		Float_t       mb_DauEPt[2];
		Float_t       mb_DauEPz[2];
		Float_t       mb_DauEpEta[2];
		Float_t       mb_DauEpPhi[2];

		Int_t         mb_DauEQ[2];
		Float_t       mb_DauEDcaXY[2];
		Float_t       mb_DauEDcaZ[2];

		Float_t       mb_DauEDca[2];
		Float_t       mb_DauEnSigE[2];
		Float_t       mb_DauELocalY[2];
		Float_t       mb_DauEInvBeta[2];

		Float_t       mb_DauEe0[2];
		Float_t       mb_DauEe[2];
		Float_t       mb_DauEpve[2];
		Int_t         mb_DauEnEta[2];
		Int_t         mb_DauEnPhi[2];
		Float_t       mb_DauEzDist[2];
		Float_t       mb_DauEphiDist[2];



		//-----------cuts--------
		Float_t     mVzCut[2];
		Float_t     mVxCut[2];
		Float_t     mVyCut[2];		
		Float_t     mVzDiffCut[2];

		Int_t       mnHitsFitCut[2];
		Int_t       mnHitsDedxCut[2];
		Float_t     mRatioCut[2];

		Float_t     mEPtCut[2];
		Float_t     mEPCut[2];
		Float_t     mEEtaCut[2];
		Float_t     mEDcaCut[2];
		Float_t     mEInvBetaCut[2];
		Float_t     mELocalYCut[2];
		Float_t     mELocalZCut[2];
		Float_t     mEnSigECut[2];

		Float_t     mPairMassCut[2];
		Float_t     mPEMassCut[2];
		Float_t     mEEPairYCut[2];

		ClassDef(StCutsTreeMaker, 1)
};

#endif
