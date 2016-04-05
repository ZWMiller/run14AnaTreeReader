#ifndef StHadronTrack_hh
#define StHadronTrack_hh

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

// Macro to control EMC variables
#define EMCON 1
class StPicoEvent;

class StHadronTrack : public TObject {
 public:
  StHadronTrack();
  ~StHadronTrack();
  StHadronTrack(StPicoDst *picoDst, StPicoTrack *t);
  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t   id() const             { return (Int_t)mId; }
  Float_t gPt() const;
  StThreeVectorF gMom() const    { return mGMom; }
  StThreeVectorF pMom() const    { return mPMom; }
  Short_t charge() const         { return (mNHitsFit>0) ? +1 : -1; }
  Int_t   nHitsFit() const       { return (mNHitsFit>0) ? (Int_t)mNHitsFit : (Int_t)(-1*mNHitsFit); }
  Int_t   nHitsDedx() const      { return (Int_t)mNHitsDedx; }
  //Int_t   nHitsMax() const       { return (Int_t)mNHitsMax; }
  Float_t nSigmaPion() const     { return (Float_t)mNSigmaPion/100.; }
  Float_t nSigmaKaon() const     { return (Float_t)mNSigmaKaon/100.; }
  //Float_t nSigmaElectron() const { return (Float_t)mNSigmaElectron/100.; }
  Float_t dca() const            { return (Float_t)mDca/10000.; }
  Float_t dcaXY() const            { return (Float_t)mDcaXY/10000.; }
  Float_t dcaZ() const            { return (Float_t)mDcaZ/10000.; }
  Bool_t isHFTTrack() const      { return mIsHftTrack; }
  //Float_t tofMatchFlag() const   { return mTofMatchFlag; }
  Float_t beta() const           { return (Float_t)mBeta/20000.; }
  Float_t localY() const         {return (Float_t)mLocalY/1000.;}
  //Float_t localZ() const          {return (Float_t)mLocalZ/1000.;}

 protected:
  UShort_t mId;               // track Id
  Short_t  mDca;              // dca * 10000
  Short_t  mDcaXY;              // dcaXY * 10000
  Short_t  mDcaZ;              // dcaZ * 10000
  StThreeVectorF mPMom;
  StThreeVectorF mGMom;
  Char_t   mNHitsFit;         // q*nHitsFit
  //Char_t   mNHitsMax;         // nHitsMax - TPC
  UChar_t  mNHitsDedx;        // nHitsDedx
  Char_t   mIsHftTrack;       // 
  Short_t  mNSigmaPion;       // nsigmaPi * 100
  Short_t  mNSigmaKaon;       // nsigmaKaon * 100
  //Short_t  mNSigmaElectron;   // nsigmaE * 100
  
  // pidTraits
  //Char_t   mTofMatchFlag;
  UShort_t mBeta;  // *20000
  Short_t  mLocalY; // *1000
  //Short_t  mLocalZ; // *1000

  friend class StPicoDst;

  ClassDef(StHadronTrack, 1)
};
#endif
