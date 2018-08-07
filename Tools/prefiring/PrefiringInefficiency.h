#ifndef PREFIRING_H
#define PREFIRING_H

/**
 * Apply MC event weight correction for trigger prefiring, parametrized by jet pt, eta
 *
 ************************************************************/

#include <string>
#include "TEfficiency.h"
#include <vector>
#include "TMath.h"
#include "Math/LorentzVector.h"
#include "TFile.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class PrefiringInefficiency {
 public:
  PrefiringInefficiency(const std::string &leptonicEfficiencyFile, const std::string &nonLeptonicEfficiencyFile);
  ~PrefiringInefficiency() {}
  
  double GetEventEfficiencyCorrection(const std::vector<LorentzVector>& jet_p4s, const bool leptonic);
  double GetEventEfficiencyCorrection(const std::vector<float>& jet_pts, const std::vector<float>& jet_etas, const bool leptonic);
  
 private:
  TEfficiency* leptonicEfficiencies;
  TEfficiency* nonLeptonicEfficiencies;
  
};

#endif
