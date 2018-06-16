#ifndef UTILS_H
#define UTILS_H
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

float DeltaR(float eta1, float eta2, float phi1, float phi2);
float DeltaR(LorentzVector v1, LorentzVector v2); 
float DeltaPhi(float phi1, float phi2);
float MT(float pt1, float phi1, float pt2, float phi2);

namespace utils {
  bool isCloseObject(const LorentzVector p1, const LorentzVector p2, const float conesize = 0.3);
  void fastOverlapRemoval(std::vector<LorentzVector>& p4js, const std::vector<LorentzVector>& p4ls, const float conesize = 0.1);
}

#endif
