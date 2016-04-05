#ifndef StPartElectronTrack_hh
#define StPartElectronTrack_hh

#include <cmath>

class StPicoTrack;
class StPicoDst;
class StPicoBTofPidTraits;
class StPicoEmcPidTraits;
class StDcaGeometry;

#include "TObject.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include <stdio.h>
#include <math.h>
#include "StEvent/StDcaGeometry.h"

class StPartElectronTrack : public TObject {
 public:
  StPartElectronTrack();
  ~StPartElectronTrack();
  StPartElectronTrack(StPicoDst *picoDst, StPicoTrack *t);
  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t   id() const             { return (Int_t)mId; }
  Float_t gPt() const;
  Float_t gEta() const;
  Float_t gPhi() const;
  StThreeVectorF gMom() const;
  StThreeVectorF const & pMom() const    { return mPMom; }
  Short_t charge() const         { return (mNHitsFit>0) ? +1 : -1; }
  Int_t   nHitsFit() const       { return (mNHitsFit>0) ? (Int_t)mNHitsFit : (Int_t)(-1*mNHitsFit); }
  //Int_t   nHitsMax() const       { return (Int_t)mNHitsMax; }
  Int_t   nHitsDedx() const      { return (Int_t)mNHitsDedx; }
  Float_t beta() const           { return (Float_t)mBeta/20000.; }
  Float_t localY() const          {return (Float_t)mLocalY/1000.;}
  //Float_t localZ() const          {return (Float_t)mLocalZ/1000.;}
  Float_t nSigmaElectron() const { return (Float_t)mNSigmaElectron/100.; }
  Float_t dca() const           { return (Float_t)mDca/10000.; }
  Float_t dcaXY() const           { return (Float_t)mDcaXY/10000.; }
  Float_t dcaZ() const           { return (Float_t)mDcaZ/10000.; }

  Bool_t isHFTTrack() const { return mIsHft; }

 protected:
  UShort_t mId;               // track Id
  StThreeVectorF mPMom;       // primary momentum, (0.,0.,0.) if none
  StThreeVectorF mGMom;       // global momentum
  Short_t  mDca;              // dca * 10000
  Short_t  mDcaXY;              // dcaXY * 10000
  Short_t  mDcaZ;              // dcaZ * 10000
  Char_t   mNHitsFit;         // q*nHitsFit
  //Char_t   mNHitsMax;         // nHitsMax
  UChar_t  mNHitsDedx;        // nHitsDedx
  Short_t  mNSigmaElectron;   // nsigmaE * 100
  Char_t   mIsHft;            // isHFT
  
  // pidTraits
  Short_t  mBeta;   // *20000
  Short_t  mLocalY; // *1000
  //Short_t  mLocalZ; // *1000

  friend class StPicoDst;

  ClassDef(StPartElectronTrack, 1)
};
inline Float_t StPartElectronTrack::gPt() const
{
  //return 1./fabs(mPar[3]);
  return mGMom.perp();
}
inline Float_t StPartElectronTrack::gEta() const
{

  //float ptt = gPt();
  //StThreeVectorF gmom(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  //return gmom.pseudoRapidity();
  return mGMom.pseudoRapidity();
}
inline Float_t StPartElectronTrack::gPhi() const
{
  //float ptt = gPt();
  //StThreeVectorF gmom(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  //return gmom.phi();
  return mGMom.phi();;
}

inline StThreeVectorF StPartElectronTrack::gMom() const
{
  //float ptt = gPt();
  //return StThreeVectorF(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  return mGMom;
}

       
#endif
