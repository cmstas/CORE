#ifndef MUONSELECTIONS_H
#define MUONSELECTIONS_H
#include "CMS3.h"
#include "TString.h"
#include "Base.h"
#include "Config.h"

//POG IDs
bool isMediumMuonPOG_forICHEP( unsigned int muIdx );
bool isLooseMuonPOG(unsigned int muIdx);
bool isMediumMuonPOG(unsigned int muIdx);
bool isTightMuonPOG(unsigned int muIdx);
bool isHighPtMuonPOG(unsigned int muIdx);
bool isBadGlobalMuon(unsigned int muIdx, bool selectClones = false);

//Main Muon ID function
bool muonID(unsigned int muIdx, id_level_t id_level);

//Tightest ID passed by muon
int muTightID(unsigned int muIdx, analysis_t analysis, int version = 1);

// tight charge requirement
int tightChargeMuon(unsigned int muIdx);

bool PassSoftMuonCut(unsigned int muIdx);

struct muIDcache {
public:
  void  setCacheValues(int idx, float miniiso, float ptratio, float ptrel) {
    idx_ = idx;
    miniiso_ = miniiso;
    ptratio_ = ptratio;
    ptrel_ = ptrel;
  }
  //make sure it was set for this electron before returning
  float getMiniiso(int idx) {assert(idx==idx_); return miniiso_;}
  float getPtratio(int idx) {assert(idx==idx_); return ptratio_;}
  float getPtrel(int idx) {assert(idx==idx_); return ptrel_;}
private:
  int idx_;
  float miniiso_;
  float ptratio_;
  float ptrel_;
};

namespace muID {
  void setCache(int idx, float miniiso, float ptratio, float ptrel);
  void unsetCache();

  // To keep updated with: https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/DataFormats/MuonReco/interface/Muon.h#L188-L212
  enum Selector {
    CutBasedIdLoose        = 1UL<< 0,  
    CutBasedIdMedium       = 1UL<< 1,  
    CutBasedIdMediumPrompt = 1UL<< 2,  // medium with IP cuts
    CutBasedIdTight        = 1UL<< 3,  
    CutBasedIdGlobalHighPt = 1UL<< 4,  // high pt muon for Z',W' (better momentum resolution)
    CutBasedIdTrkHighPt    = 1UL<< 5,  // high pt muon for boosted Z (better efficiency)
    PFIsoVeryLoose         = 1UL<< 6,  // reliso<0.40
    PFIsoLoose             = 1UL<< 7,  // reliso<0.25
    PFIsoMedium            = 1UL<< 8,  // reliso<0.20
    PFIsoTight             = 1UL<< 9,  // reliso<0.15
    PFIsoVeryTight         = 1UL<<10,  // reliso<0.10
    TkIsoLoose             = 1UL<<11,  // reliso<0.10
    TkIsoTight             = 1UL<<12,  // reliso<0.05
    SoftCutBasedId         = 1UL<<13,  
    SoftMvaId              = 1UL<<14,  
    MvaLoose               = 1UL<<15,  
    MvaMedium              = 1UL<<16,  
    MvaTight               = 1UL<<17,
    MiniIsoLoose           = 1UL<<18,  // reliso<0.40
    MiniIsoMedium          = 1UL<<19,  // reliso<0.20
    MiniIsoTight           = 1UL<<20,  // reliso<0.10
    MiniIsoVeryTight       = 1UL<<21,  // reliso<0.05
    TriggerIdLoose         = 1UL<<22,  // robust selector for HLT
    InTimeMuon             = 1UL<<23,
    PFIsoVeryVeryTight     = 1UL<<24,  // reliso<0.05
    MultiIsoLoose          = 1UL<<25,  // miniIso with ptRatio and ptRel
    MultiIsoMedium         = 1UL<<26   // miniIso with ptRatio and ptRel
  };
}

// Get muon identification info from MiniAOD (availability start from 94X)
bool passesMuonPOG(muID::Selector id_type, int idx);

// check CMS3 version to see which c++ type is stored in the ntuples for mus_gfit_ndof
int get_mus_gfit_ndof( unsigned int muIdx );

#endif
