#include <vector>
#include "utils.h"
#include "TMath.h"

using namespace std;
const float PI = TMath::Pi();

float DeltaR(float eta1, float eta2, float phi1, float phi2){
  float dEta = eta1 - eta2;
  float dPhi = DeltaPhi(phi1, phi2);
  return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

float DeltaR(LorentzVector v1, LorentzVector v2){
  return DeltaR(v1.eta(), v2.eta(), v1.phi(), v2.phi()); 
}

float DeltaPhi(float phi1, float phi2){
  float dPhi = phi1 - phi2;
  while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
  return fabs(dPhi);
}

float MT(float pt1, float phi1, float pt2, float phi2){
  return sqrt( 2 * pt1 * pt2 * ( 1 - cos( phi1 - phi2 ) ) );
}

bool utils::isCloseObject(const LorentzVector p1, const LorentzVector p2, const float conesize) {
  float deltaEta = fabs(p1.eta() - p2.eta());
  if (deltaEta > conesize) return false;
  float deltaPhi = fabs(p1.phi() - p2.phi());
  if (deltaPhi > PI) deltaPhi = 2*PI - deltaPhi;
  if (deltaPhi > conesize) return false;
  float deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
  if (deltaR2 > conesize*conesize) return false;

  return true;
}

void utils::fastOverlapRemoval(vector<LorentzVector>& p4jets, const vector<LorentzVector>& p4leps, const float conesize) {
  for (auto lep : p4leps) {
    float minDR = 999;
    auto rm_it = p4jets.end();
    for (auto ijet = p4jets.begin(); ijet < p4jets.end(); ++ijet) {
      auto& jet = *ijet;
      float deltaEta = fabs(jet.eta() - lep.eta());
      if (deltaEta > conesize) continue;
      float deltaPhi = fabs(jet.phi() - lep.phi());
      if (deltaPhi > PI) deltaPhi = 2*PI - deltaPhi;
      if (deltaPhi > conesize) continue;
      float deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
      if (deltaR > conesize) continue;
      if (deltaR < minDR) {
        minDR = deltaR;
        rm_it = ijet;
      }
    }
    if (minDR < conesize)
      p4jets.erase(rm_it);
  }
}
