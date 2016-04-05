#include "StElectronTrack.h"
#include "StMessMgr.h"
#include "TVector2.h"
#include "TMath.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StBTofUtil/tofPathLength.hh"
#include "PhysicalConstants.h"
#include <climits>

ClassImp(StElectronTrack)

//----------------------------------------------------------------------------------
StElectronTrack::StElectronTrack() : mId(0), mPMom(0., 0., 0.), mGMom(0.,0.,0.), mDedx(0),
   mNHitsFit(0), mNHitsDedx(0), 
   mNSigmaElectron(32768), mMap0(0),
   mCurv(0), mDip(0), mPhase(0), mOrigin(0,0,0), mH(0), 
   mBeta(0), mLocalY(32768),
   mBTOWADC0(0), mBTOWE0(0), mBTOWE(0),
   mBEMCDistZ(32768), mBEMCDistPhi(32768), mBSMDNEta(0), mBSMDNPhi(0),
   mBTOWId(0),  mEmcTrgId(-1)
{

}

//----------------------------------------------------------------------------------
StElectronTrack::StElectronTrack(StPicoDst *picoDst, StPicoTrack* t)
{
   mId        = (UShort_t)t->id();
   //mChi2      = (t->chi2() * 1000. > 65536) ? 65536 : (UShort_t)(TMath::Nint(t->chi2() * 1000.));
   mPMom = t->pMom();
   int q      = t->charge();
   mDedx      = (t->dEdx()*1000. > 65536) ? 65536 : (UShort_t)(TMath::Nint(t->dEdx()*1000.));
   mNHitsFit  = t->nHitsFit()*q;
   //mNHitsMax  = t->nHitsMax();
   mNHitsDedx = (UChar_t)(t->nHitsDedx());
   //mNSigmaPion = (fabs(t->nSigmaPion() * 100.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaPion() * 100.));
   //mNSigmaKaon = (fabs(t->nSigmaKaon() * 100.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaKaon() * 100.));
   //mNSigmaProton = (fabs(t->nSigmaProton() * 100.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaProton() * 100.));
   mNSigmaElectron = (fabs(t->nSigmaElectron() * 100.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaElectron() * 100.));

   mMap0 = t->map0(); // see hitMap definition in StTrackTopologyMap
   //mMap1 = t->map1();
   //const float* params = t->params();
   //const float* errMatrix = t->errMatrix();
   //for (int i = 0; i < 6; i++) mPar[i] = params[i];
   //for (int i = 0; i < 15; i++) mErrMatrix[i] = errMatrix[i];

   StThreeVectorF vertexPos = picoDst->event()->primaryVertex();
   StPhysicalHelixD helix = t->helix();
   mGMom = t->gMom(vertexPos,picoDst->event()->bField());
   mOrigin = helix.origin();
   mCurv = helix.curvature();
   mDip = helix.dipAngle();
   mPhase = helix.phase();
   mH = helix.h();
   StThreeVectorF dcaPoint = helix.at(helix.pathLength(vertexPos.x(), vertexPos.y()));
   float dcaZ = (dcaPoint.z() - vertexPos.z())*10000.;
   float dcaXY = (helix.geometricSignedDistance(vertexPos.x(),vertexPos.y()))*10000.;
   mDcaZ = dcaZ>32768?32768:(Short_t)dcaZ;
   mDcaXY = dcaXY>32768?32768:(Short_t)dcaXY;

   double thePath = helix.pathLength(vertexPos);
   StThreeVectorF dcaPos = helix.at(thePath);
   mDca = fabs((dcaPos-vertexPos).mag()*10000.)>32768? 32768: (Short_t)((dcaPos-vertexPos).mag()*10000.);

   int index2TofPid = t->bTofPidTraitsIndex();
   if (index2TofPid>=0){
      StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2TofPid);
      //mTofMatchFlag = tofPid->btofMatchFlag();
      Float_t mom = mPMom.mag();
      mLocalY = tofPid->btofYLocal()*1000;
      //mLocalZ = tofPid->btofZLocal()*1000;
      Float_t beta = tofPid->btofBeta();
      if(beta<1e-4||beta>=(USHRT_MAX-1)/20000){
         Float_t tof = tofPid->btof();
         StThreeVectorF btofHitPos = tofPid->btofHitPos();
         float L = tofPathLength(&vertexPos, &btofHitPos, helix.curvature()); 
         if(tof>0) beta = L/(tof*(c_light/1.0e9));
      }
      mBeta = (UShort_t)(beta*20000);
   }

   int index2EmcPid = t->emcPidTraitsIndex();
   if (index2EmcPid>=0){
      StPicoEmcPidTraits *emcPid = picoDst->emcPidTraits(index2EmcPid);
      mBTOWADC0 = emcPid->adc0();
      mBTOWE0 = emcPid->e0()*1000;
      mBTOWE = emcPid->e()*1000;
      mBEMCDistZ = emcPid->zDist()*1000;
      mBEMCDistPhi = emcPid->phiDist()*10000;
      mBSMDNEta = emcPid->nEta();
      mBSMDNPhi = emcPid->nPhi();

      mBTOWId = emcPid->btowId();
      mEmcTrgId = -1;
   }else{
      mBTOWADC0 = 0;
      mBTOWE0 = 0;
      mBTOWE = 0;
      mBEMCDistZ = 32768;
      mBEMCDistPhi = 32768;
      mBSMDNEta = 0;
      mBSMDNPhi = 0;

      mBTOWId = 0;
      mEmcTrgId = -1;

   }
}

//----------------------------------------------------------------------------------
StElectronTrack::~StElectronTrack()
{
   /* noop */
}
//----------------------------------------------------------------------------------
void StElectronTrack::Print(const Char_t *option) const
{
   if (strcmp(option, "tpc") == 0 || strcmp(option, "") == 0)
   {
      LOG_INFO << "id=" << id() <<" charge = "<<charge()
         << endm;
      LOG_INFO << "pMom=" << pMom() << endm;
      LOG_INFO << "gpt =" << gPt() << " gEta = "<<gEta()<<" gPhi = "<<gPhi()<<endm;
      LOG_INFO << " nHitsFit = " << nHitsFit() << " nHitsdEdx = " << nHitsDedx() << endm;
      LOG_INFO << " nSigma E = " << nSigmaElectron() << endm;
      LOG_INFO << " beta = "<<beta()<<" dca = "<<dca()<<endm;
   }
}
