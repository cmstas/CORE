#include "PrefiringInefficiency.h"
#include <iostream>
#include <exception>

PrefiringInefficiency::PrefiringInefficiency(const std::string &leptonicEfficiencyFile, const std::string &nonLeptonicEfficiencyFile) {
  TFile* leptonicFile = new TFile(leptonicEfficiencyFile.c_str(), "READ");
  if (leptonicFile->IsZombie()) {
    std::cout << "[Trigger Prefiring] File " << leptonicEfficiencyFile << " not found." << std::endl;
    throw std::exception();
  }
  leptonicEfficiencies = (TEfficiency*) leptonicFile->Get("prefireEfficiencyMap")->Clone("leptonicEfficiencies");
  leptonicFile->Close();

  TFile* nonLeptonicFile = new TFile(nonLeptonicEfficiencyFile.c_str(), "READ");
  if (nonLeptonicFile->IsZombie()) {
    std::cout << "[Trigger Prefiring] File " << leptonicEfficiencyFile << " not found." << std::endl;
    throw std::exception();
  }
  nonLeptonicEfficiencies = (TEfficiency*) nonLeptonicFile->Get("prefireEfficiencyMap")->Clone("nonLeptonicEfficiencies");  
  nonLeptonicFile->Close();
}

double PrefiringInefficiency::GetEventEfficiencyCorrection(const std::vector<LorentzVector>& jet_p4s, const bool leptonic) {
  TEfficiency* efficiencies = leptonic ? leptonicEfficiencies : nonLeptonicEfficiencies;
  double weight = 1.0;
  for (std::vector<LorentzVector>::const_iterator jet_it = jet_p4s.begin(); jet_it != jet_p4s.end(); jet_it++) {
    double efficiency = efficiencies->GetEfficiency(efficiencies->FindFixBin(fabs(jet_it->eta()),jet_it->pt()));
    weight *= (1-efficiency);
  }
  return weight;
}

double PrefiringInefficiency::GetEventEfficiencyCorrection(const std::vector<float>& jet_pts, const std::vector<float>& jet_etas, const bool leptonic) {
  if (jet_pts.size() != jet_etas.size()) {
    std::cout << "[Triggering Prefiring] Jet pt and eta vector must be the same size." << std::endl;
    throw std::exception();
  }
  TEfficiency* efficiencies = leptonic ? leptonicEfficiencies : nonLeptonicEfficiencies;
  double weight = 1.0;
  for (std::vector<float>::const_iterator jet_pt_it = jet_pts.begin(); jet_pt_it != jet_pts.end(); jet_pt_it++) {
    const float jet_eta = fabs(jet_etas.at(jet_pt_it - jet_pts.begin()));
    double efficiency = efficiencies->GetEfficiency(efficiencies->FindFixBin(jet_eta, *jet_pt_it));
    weight *= (1-efficiency);
  }
  return weight;
}
