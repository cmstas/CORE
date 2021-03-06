#include "Math/VectorUtil.h"
#include "IsolationTools.h"
#include "Tools/JetCorrector.h"

using namespace std;
using namespace tas;

bool passMultiIsoCuts(float cutMiniIso, float cutPtRatio, float cutPtRel, float miniIsoValue, float ptRatioValue, float ptRelValue){
  return (miniIsoValue < cutMiniIso && (ptRatioValue>cutPtRatio || ptRelValue > cutPtRel));
}

bool passMultiIso(int id, int idx, float cutMiniIso, float cutPtRatio, float cutPtRel, int eaversion, int whichCorr){
  const LorentzVector& lep_p4 = abs(id)==11 ? els_p4().at(idx) : mus_p4().at(idx);
  const LorentzVector& jet_p4 = closestJet(lep_p4,0.4,3.0,whichCorr);
  float miniIso = abs(id)==11 ? elMiniRelIsoCMS3_EA(idx,eaversion) : muMiniRelIsoCMS3_EA(idx,eaversion);
  float closeJetPt = jet_p4.pt();
  float ptratio = ( closeJetPt>0. ? lep_p4.pt()/closeJetPt : 1.);
  float ptrel = ptRel(lep_p4, jet_p4, true);
  return passMultiIsoCuts(cutMiniIso, cutPtRatio, cutPtRel, miniIso, ptratio, ptrel);
}

float muMiniRelIsoCMS3_charged(unsigned int idx) {
  float absiso = mus_miniIso_ch().at(idx);
  return absiso/(mus_p4().at(idx).pt());
}

bool passPtRel(int id, int idx, float cut, bool subtractLep, int whichCorr) {
  return getPtRel(id, idx, subtractLep, whichCorr) > cut;
}

float getPtRel(int id, int idx, bool subtractLep, int whichCorr) {
  LorentzVector lep_p4 = abs(id)==11 ? els_p4().at(idx) : mus_p4().at(idx);
  LorentzVector jet_p4 = closestJet(lep_p4,0.4,3.0,whichCorr);
  return ptRel(lep_p4, jet_p4, subtractLep);
}

float ptRel(const LorentzVector& lepp4, const LorentzVector& jetp4, bool subtractLep) {
  if (jetp4.pt()<=0.) return 0.;
  LorentzVector jp4 = jetp4;
  if (subtractLep) jp4-=lepp4;
  float dot = lepp4.Vect().Dot( jp4.Vect() );
  float ptrel = lepp4.P2() - dot*dot/jp4.P2();
  ptrel = ptrel>0 ? sqrt(ptrel) : 0.0;
  return ptrel;
}

int closestJetIdx(const LorentzVector& lep_p4, float dRmin, float maxAbsEta){
  int closestIdx = -1;
  for (unsigned int pfjidx=0; pfjidx<pfjets_p4().size(); ++pfjidx){
    LorentzVector jet_p4 = pfjets_p4().at(pfjidx);
    if (fabs(jet_p4.eta()) > maxAbsEta) continue;
    float tmp_dRmin = ROOT::Math::VectorUtil::DeltaR(jet_p4, lep_p4);
    if (tmp_dRmin < dRmin){
      closestIdx = pfjidx;
      dRmin = tmp_dRmin;
    }
  }
  return closestIdx;
}

LorentzVector closestJet(const LorentzVector& lep_p4, float dRmin, float maxAbsEta, int whichCorr){
  int closestIdx = closestJetIdx(lep_p4,dRmin,maxAbsEta);
  if (closestIdx < 0) return LorentzVector();
  LorentzVector jet = pfjets_p4().at(closestIdx);
  LorentzVector rawjet = pfjets_p4().at(closestIdx)*pfjets_undoJEC().at(closestIdx);

  if (whichCorr == 0) return jet;

  //Calculate JEC
  if (gconf.jet_corrector_L1 == 0){
      throw std::runtime_error("jet corrector for L1 not initialized but needed for closestJet in IsolationTools");
  }
  if (gconf.jet_corrector_L2L3 == 0){
      throw std::runtime_error("jet corrector for L2L3 not initialized but needed for closestJet in IsolationTools");
  }
  gconf.jet_corrector_L1->setJetEta(rawjet.eta());
  gconf.jet_corrector_L1->setJetPt(rawjet.pt());
  gconf.jet_corrector_L1->setJetA(tas::pfjets_area().at(closestIdx));
  gconf.jet_corrector_L1->setRho(tas::evt_fixgridfastjet_all_rho());
  float JEC1 = gconf.jet_corrector_L1->getCorrection();
  gconf.jet_corrector_L2L3->setJetEta(rawjet.eta());
  gconf.jet_corrector_L2L3->setJetPt(rawjet.pt()*JEC1);
  gconf.jet_corrector_L2L3->setJetA(tas::pfjets_area().at(closestIdx));
  gconf.jet_corrector_L2L3->setRho(tas::evt_fixgridfastjet_all_rho());
  float JEC2 = gconf.jet_corrector_L2L3->getCorrection();

  if (whichCorr == 1) return rawjet*JEC1;
  return (rawjet*JEC1 - lep_p4)*JEC2 + lep_p4;
}

float getMiniDR(float pt) {
  float dr = 0.2;
  if (pt>50) dr = 10./pt;
  if (pt>200) dr = 0.05;
  return dr;
}

float muRelIso03DB(unsigned int muIdx){
  float chiso     = mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR03_pf_PhotonEt().at(muIdx);
  float deltaBeta = mus_isoR03_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(mus_p4().at(muIdx).pt());
}

float muRelIso04DB(unsigned int muIdx){
  float chiso     = mus_isoR04_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR04_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR04_pf_PhotonEt().at(muIdx);
  float deltaBeta = mus_isoR04_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(mus_p4().at(muIdx).pt());
}

float muRelIso03(unsigned int muIdx, analysis_t analysis){
  if (analysis == WW  ) return muRelIso03EA(muIdx);
  if (analysis == SS  ) return muRelIso03EA(muIdx,gconf.ea_version);
  if (analysis == ZMET) return muRelIso03EA(muIdx);
  return muRelIso03DB(muIdx);
}

float muRelIso04(unsigned int muIdx, analysis_t analysis){
  return muRelIso04DB(muIdx);
}

float muEA03(unsigned int muIdx, int version){
  float ea = 0.;
  if (version==0) {
    //PHYS14 version
    if      (fabs(mus_p4().at(muIdx).eta())<=0.800) ea = 0.0913;
    else if (fabs(mus_p4().at(muIdx).eta())<=1.300) ea = 0.0765;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.000) ea = 0.0546;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.200) ea = 0.0728;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.500) ea = 0.1177;
  } else if (version==3 || version==4) {
      // Fall17 https://github.com/cms-data/PhysicsTools-NanoAOD/blob/master/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt
    if      (fabs(mus_p4().at(muIdx).eta())<=0.800) ea = 0.0566;
    else if (fabs(mus_p4().at(muIdx).eta())<=1.300) ea = 0.0562;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.000) ea = 0.0363;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.200) ea = 0.0119;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.500) ea = 0.0064;
  } else if(version==1 || version==2){
    //Spring15 version
    //For version 2 it is same as Spring 15. Why? I don't know. I think Dylan probably knows about this.)
    if      (fabs(mus_p4().at(muIdx).eta())<=0.800) ea = 0.0735;
    else if (fabs(mus_p4().at(muIdx).eta())<=1.300) ea = 0.0619;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.000) ea = 0.0465;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.200) ea = 0.0433;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.500) ea = 0.0577;
  }
  else{
      // invalid value
      cerr << "[CORE:IsolationTools:muEA03] WARNING! invalid EA version value. (Did you forget to change your gconf.ea_version?)" << endl;
  }

  return ea;
}

float muRelIso03EA(unsigned int muIdx, int eaversion, bool include_leptons){
  float chiso     = mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR03_pf_PhotonEt().at(muIdx);
  float ea = muEA03(muIdx, eaversion);
  if (include_leptons){
    chiso = mus_isoR03_pf_ChargedParticlePt().at(muIdx);
  }
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - evt_fixgridfastjet_all_rho() * ea);
  return absiso/(mus_p4().at(muIdx).pt());
}

float muRelIso03_noCorr(unsigned int muIdx){
  float chiso  = mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso  = mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso  = mus_isoR03_pf_PhotonEt().at(muIdx);
  float absiso = chiso + nhiso + emiso;
  return absiso/(mus_p4().at(muIdx).pt());
}

float muRelIsoCustomCone(unsigned int muIdx, float dr, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor, float mindr, int eaversion, bool include_leptons){
  cout << "WARNING: shouldn't use muRelIsoCustomCone in CMS4! It depends on entire set of PFCands." << endl;

  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  float correction = 0.;
  float deltaBeta = 0.;
  float deadcone_ch = 0.0001;
  float deadcone_pu = 0.01;
  float deadcone_ph = 0.01;
  float deadcone_nh = 0.01;

  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    float thisDR = fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(i),mus_p4().at(muIdx)));
    if (thisDR<mindr) continue;
    if (thisDR>dr) continue;
    if ( fabs(pfcands_particleId().at(i))==211 || (include_leptons && ((fabs(pfcands_particleId().at(i)) == 11) || (fabs(pfcands_particleId().at(i)) == 13 && dr > 0.0005))) ) {
      if (pfcands_fromPV().at(i) > 1 && (!useVetoCones || thisDR > deadcone_ch) ) chiso+=pfcands_p4().at(i).pt();
      else if (useDBcor && pfcands_fromPV().at(i) <= 1 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_pu)) deltaBeta+=pfcands_p4().at(i).pt();
    }
    if (fabs(pfcands_particleId().at(i))==130 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_nh)) nhiso+=pfcands_p4().at(i).pt();
    if (fabs(pfcands_particleId().at(i))==22 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_ph)) emiso+=pfcands_p4().at(i).pt();
  }

  if (useDBcor) correction = 0.5 * deltaBeta;
  else if (useEAcor) correction = evt_fixgridfastjet_all_rho() * muEA03(muIdx,eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - correction);
  return absiso/(mus_p4().at(muIdx).pt());
}

[[deprecated("This function may not be correct in cms4 (needs all PFCands)")]]
float muMiniRelIso(unsigned int idx, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor){
  float pt = mus_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return  muRelIsoCustomCone(idx,dr,useVetoCones,ptthresh,useDBcor,useEAcor);
}

float muMiniRelIsoCMS3_DB(unsigned int idx) {
  float correction = 0.5 * mus_miniIso_db().at(idx);
  float absiso = mus_miniIso_ch().at(idx) + std::max(float(0.0), mus_miniIso_nh().at(idx) + mus_miniIso_em().at(idx) - correction);
  return absiso/(mus_p4().at(idx).pt());
}

float muMiniRelIsoCMS3_EA(unsigned int idx, int eaversion) {
  float pt = mus_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  float rho = evt_fixgridfastjet_all_rho();
  // if (eaversion!=0) rho = evt_fixgridfastjet_centralneutral_rho();
  float correction = rho * muEA03(idx,eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = mus_miniIso_ch().at(idx) + std::max(float(0.0), mus_miniIso_nh().at(idx) + mus_miniIso_em().at(idx) - correction);
  return absiso/(mus_p4().at(idx).pt());
}

[[deprecated("This function may not be correct in cms4 (needs all PFCands)")]]
float muRelIsoAn04(unsigned int idx, bool useDBcor) {
  float pt = mus_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return muRelIsoCustomCone(idx, 0.4, false, 0.0, useDBcor, false, dr);
}

float muPtRatio(unsigned int idx)
{
  const LorentzVector& temp_jet_p4 = closestJet(cms3.mus_p4()[idx], 0.4, 3.0, 2);
  float closeJetPt = temp_jet_p4.pt();
  float ptratio = (closeJetPt > 0. ? cms3.mus_p4()[idx].pt() / closeJetPt : 1.);
  return ptratio;
}

float eleRelIso03(unsigned int elIdx, analysis_t analysis){
  if (analysis == HAD ) return eleRelIso03DB(elIdx);
  if (analysis == STOP) return eleRelIso03DB(elIdx);
  if (analysis == WW  ) return eleRelIso03_90ContEA(elIdx);
  if (analysis == SS  ) return eleRelIso03EA(elIdx,gconf.ea_version);
  if (analysis == ZMET) return eleRelIso03EA(elIdx);

  else return eleRelIso03EA(elIdx,0);
}

float eleRelIso03DB(unsigned int elIdx){
  float chiso     = els_pfChargedHadronIso().at(elIdx);
  float nhiso     = els_pfNeutralHadronIso().at(elIdx);
  float emiso     = els_pfPhotonIso().at(elIdx);
  float deltaBeta = els_pfPUIso().at(elIdx);
  float absiso    = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(els_p4().at(elIdx).pt());
}

float elEA03(unsigned int elIdx, int version) {
  float ea = 0.;
  if (version==0) {
    //PHYS14 version
    if      (fabs(els_p4().at(elIdx).eta())<=0.800) ea = 0.1013;
    else if (fabs(els_p4().at(elIdx).eta())<=1.300) ea = 0.0988;
    else if (fabs(els_p4().at(elIdx).eta())<=2.000) ea = 0.0572;
    else if (fabs(els_p4().at(elIdx).eta())<=2.200) ea = 0.0842;
    else if (fabs(els_p4().at(elIdx).eta())<=2.500) ea = 0.1530;
  }
  else if (version==1){
    //Spring15 version
    if      (fabs(els_etaSC().at(elIdx))<=1.000) ea = 0.1752;
    else if (fabs(els_etaSC().at(elIdx))<=1.479) ea = 0.1862;
    else if (fabs(els_etaSC().at(elIdx))<=2.000) ea = 0.1411;
    else if (fabs(els_etaSC().at(elIdx))<=2.200) ea = 0.1534;
    else if (fabs(els_etaSC().at(elIdx))<=2.300) ea = 0.1903;
    else if (fabs(els_etaSC().at(elIdx))<=2.400) ea = 0.2243;
    else if (fabs(els_etaSC().at(elIdx))<=2.500) ea = 0.2687;
  }
  else if (version==2){
    //Spring16 version
    if      (fabs(els_etaSC().at(elIdx))<=1.000) ea = 0.1703;
    else if (fabs(els_etaSC().at(elIdx))<=1.479) ea = 0.1715;
    else if (fabs(els_etaSC().at(elIdx))<=2.000) ea = 0.1213;
    else if (fabs(els_etaSC().at(elIdx))<=2.200) ea = 0.1230;
    else if (fabs(els_etaSC().at(elIdx))<=2.300) ea = 0.1635;
    else if (fabs(els_etaSC().at(elIdx))<=2.400) ea = 0.1937;
    else if (fabs(els_etaSC().at(elIdx))<=2.500) ea = 0.2393;
  }
  else if (version==3) {
    //Fall17 92X version https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt
    if      (fabs(els_etaSC().at(elIdx))<=1.000) ea = 0.1566;
    else if (fabs(els_etaSC().at(elIdx))<=1.479) ea = 0.1626;
    else if (fabs(els_etaSC().at(elIdx))<=2.000) ea = 0.1073;
    else if (fabs(els_etaSC().at(elIdx))<=2.200) ea = 0.0854;
    else if (fabs(els_etaSC().at(elIdx))<=2.300) ea = 0.1051;
    else if (fabs(els_etaSC().at(elIdx))<=2.400) ea = 0.1204;
    else if (fabs(els_etaSC().at(elIdx))<=2.500) ea = 0.1524;
  }
  else if (version==4) {
    //Fall17 94X version https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
    if      (fabs(els_etaSC().at(elIdx))<=1.000) ea = 0.1440;
    else if (fabs(els_etaSC().at(elIdx))<=1.479) ea = 0.1562;
    else if (fabs(els_etaSC().at(elIdx))<=2.000) ea = 0.1032;
    else if (fabs(els_etaSC().at(elIdx))<=2.200) ea = 0.0859;
    else if (fabs(els_etaSC().at(elIdx))<=2.300) ea = 0.1116;
    else if (fabs(els_etaSC().at(elIdx))<=2.400) ea = 0.1321;
    else if (fabs(els_etaSC().at(elIdx))<=2.500) ea = 0.1654;
  }
  else {
      // invalid value
      cerr << "[CORE:IsolationTools:elEA03] WARNING! invalid EA version value. (Did you forget to change your gconf.ea_version?)" << endl;
  }
  return ea;
}

float el90ContEA03(unsigned int elIdx) {
  float ea = 0.;
  if      (fabs(els_p4().at(elIdx).eta())<=1.000) ea = 0.1752;
  else if (fabs(els_p4().at(elIdx).eta())<=1.479) ea = 0.1862;
  else if (fabs(els_p4().at(elIdx).eta())<=2.000) ea = 0.1411;
  else if (fabs(els_p4().at(elIdx).eta())<=2.200) ea = 0.1534;
  else if (fabs(els_p4().at(elIdx).eta())<=2.300) ea = 0.1903;
  else if (fabs(els_p4().at(elIdx).eta())<=2.400) ea = 0.2243;
  else if (fabs(els_p4().at(elIdx).eta())<=2.500) ea = 0.2687;
  return ea;
}

float eleRelIso03EA(unsigned int elIdx, int eaversion, bool include_leptons){
  float chiso = els_pfChargedHadronIso().at(elIdx);
  float nhiso = els_pfNeutralHadronIso().at(elIdx);
  float emiso = els_pfPhotonIso().at(elIdx);
  float ea    = elEA03(elIdx, eaversion);
  float lepiso = 0;
  if (include_leptons){
    for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
      float thisDR = fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(i),els_p4().at(elIdx)));
      if (thisDR>0.3 || thisDR<0.0005) // To avoid adding itself
        continue;
      if (fabs(pfcands_particleId().at(i))==11 && pfcands_fromPV().at(i) > 1){
        lepiso += pfcands_p4().at(i).pt();
      }
      if (fabs(pfcands_particleId().at(i))==13 && pfcands_fromPV().at(i) > 1){
        lepiso += pfcands_p4().at(i).pt();
      }
    }
  }
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - evt_fixgridfastjet_all_rho() * ea) + lepiso;
  return absiso/(els_p4().at(elIdx).pt());
}

float eleRelIso03_90ContEA(unsigned int elIdx){
  float chiso = els_pfChargedHadronIso().at(elIdx);
  float nhiso = els_pfNeutralHadronIso().at(elIdx);
  float emiso = els_pfPhotonIso().at(elIdx);
  float ea    = el90ContEA03(elIdx);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - evt_fixgridfastjet_all_rho() * ea);
  return absiso/(els_p4().at(elIdx).pt());
}

float eleRelIso03_noCorr(unsigned int elIdx){
  float chiso = els_pfChargedHadronIso().at(elIdx);
  float nhiso = els_pfNeutralHadronIso().at(elIdx);
  float emiso = els_pfPhotonIso().at(elIdx);
  float absiso = chiso + nhiso + emiso;
  return absiso/(els_p4().at(elIdx).pt());
}

float elRelIsoCustomCone(unsigned int elIdx, float dr, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor, float mindr, int eaversion, bool include_leptons){
  cout << "WARNING: shouldn't use elRelIsoCustomCone in CMS4! It depends on entire set of PFCands." << endl;

  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  float correction = 0.;
  float deltaBeta = 0.;
  float deadcone_ch = 0.;
  float deadcone_pu = 0.;
  float deadcone_ph = 0.;
  // veto cones only in the endcap for electrons
  if (useVetoCones && fabs(els_etaSC().at(elIdx)) > 1.479) {
    deadcone_ch = 0.015;
    deadcone_pu = 0.015;
    deadcone_ph = 0.08;
  }
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    float thisDR = fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(i),els_p4().at(elIdx)));
    if ( thisDR<mindr ) continue;
    if ( thisDR>dr ) continue;
    if ( fabs(pfcands_particleId().at(i))==211  || (include_leptons && ((fabs(pfcands_particleId().at(i)) == 11 && dr > 0.0005) || (fabs(pfcands_particleId().at(i)) == 13)))) {
      if (pfcands_fromPV().at(i) > 1 && (!useVetoCones || thisDR > deadcone_ch) ) chiso+=pfcands_p4().at(i).pt();
      else if (useDBcor && pfcands_fromPV().at(i) <= 1 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_pu)) deltaBeta+=pfcands_p4().at(i).pt();
    }
    if ( fabs(pfcands_particleId().at(i))==130 && (pfcands_p4().at(i).pt() > ptthresh) ) nhiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==22 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_ph) ) emiso+=pfcands_p4().at(i).pt();
  }
  if (useDBcor) correction = 0.5 * deltaBeta;
  else if (useEAcor) correction = evt_fixgridfastjet_all_rho() * elEA03(elIdx, eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - correction);
  return absiso/(els_p4().at(elIdx).pt());
}

[[deprecated("This function may not be correct in cms4 (needs all PFCands)")]]
float elMiniRelIso(unsigned int idx, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor){

  float pt = els_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return elRelIsoCustomCone(idx,dr,useVetoCones,ptthresh,useDBcor,useEAcor);
}

float elMiniRelIsoCMS3_DB(unsigned int idx) {
  float correction = 0.5 * els_miniIso_db().at(idx);
  float absiso = els_miniIso_ch().at(idx) + std::max(float(0.0), els_miniIso_nh().at(idx) + els_miniIso_em().at(idx) - correction);
  return absiso/(els_p4().at(idx).pt());
}

float elMiniRelIsoCMS3_EA(unsigned int idx, int eaversion) {
  float pt = els_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  float rho = evt_fixgridfastjet_all_rho();
  // if (eaversion!=0) rho = evt_fixgridfastjet_centralneutral_rho();
  float correction = rho * elEA03(idx, eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = els_miniIso_ch().at(idx) + std::max(float(0.0), els_miniIso_nh().at(idx) + els_miniIso_em().at(idx) - correction);
  return absiso/(els_p4().at(idx).pt());
}

float elMiniRelIsoCMS3_charged(unsigned int idx) {
  float absiso = els_miniIso_ch().at(idx);
  return absiso/(els_p4().at(idx).pt());
}

[[deprecated("This function may not be correct in cms4 (needs all PFCands)")]]
float elRelIsoAn04(unsigned int idx, bool useDBcor) {
  float pt = els_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return elRelIsoCustomCone(idx, 0.4, false, 0.0, useDBcor, false, dr);
}

float elPtRatio(unsigned int idx)
{
  const LorentzVector& temp_jet_p4 = closestJet(cms3.els_p4()[idx], 0.4, 3.0, 2);
  float closeJetPt = temp_jet_p4.pt();
  float ptratio = (closeJetPt > 0. ? cms3.els_p4()[idx].pt() / closeJetPt : 1.);
  return ptratio;
}

//////////////////////
// Photon Isolation //
//////////////////////

float photon_CHEA03( int photonIdx, int version )
{
  float aeta = fabs(cms3.photons_p4().at(photonIdx).eta());
  float EA = -999;

  if(      version == 0 ){
    // from https://indico.cern.ch/event/369239/contribution/2/attachments/1134693/1623149/spring15_pcb.pdf
    // see also: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2Archive#Selection_implementation_details
    if(      aeta < 1.0   ) EA = 0.0157;
    else if( aeta < 1.479 ) EA = 0.0143;
    else if( aeta < 2.0   ) EA = 0.0115;
    else if( aeta < 2.2   ) EA = 0.0094;
    else if( aeta < 2.3   ) EA = 0.0095;
    else if( aeta < 2.4   ) EA = 0.0068;
    else if( aeta > 2.4   ) EA = 0.0053;
  }
  else if( version == 2 ){
    // Spring16: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
    if(      aeta < 1.0   ) EA = 0.0360;
    else if( aeta < 1.479 ) EA = 0.0377;
    else if( aeta < 2.0   ) EA = 0.0306;
    else if( aeta < 2.2   ) EA = 0.0283;
    else if( aeta < 2.3   ) EA = 0.0254;
    else if( aeta < 2.4   ) EA = 0.0217;
    else if( aeta > 2.4   ) EA = 0.0167;
  }
  else if( version == 4 ){
    // Fall17 V2: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
    if(      aeta < 1.0   ) EA = 0.0112;
    else if( aeta < 1.479 ) EA = 0.0108;
    else if( aeta < 2.0   ) EA = 0.0106;
    else if( aeta < 2.2   ) EA = 0.01002;
    else if( aeta < 2.3   ) EA = 0.0098;
    else if( aeta < 2.4   ) EA = 0.0089;
    else if( aeta > 2.4   ) EA = 0.0087;
  }

  return EA;
}

float photon_NHEA03( int photonIdx, int version )
{
  float aeta = fabs(cms3.photons_p4().at(photonIdx).eta());
  float EA = -999;

  if(      version == 0 ){
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2Archive#Selection_implementation_details
    if(      aeta < 1.0   ) EA = 0.0143;
    else if( aeta < 1.479 ) EA = 0.0210;
    else if( aeta < 2.0   ) EA = 0.0148;
    else if( aeta < 2.2   ) EA = 0.0082;
    else if( aeta < 2.3   ) EA = 0.0124;
    else if( aeta < 2.4   ) EA = 0.0186;
    else if( aeta > 2.4   ) EA = 0.0320;
  }
  else if( version == 2 ){
    // Spring16: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
    if(      aeta < 1.0   ) EA = 0.0597;
    else if( aeta < 1.479 ) EA = 0.0807;
    else if( aeta < 2.0   ) EA = 0.0629;
    else if( aeta < 2.2   ) EA = 0.0197;
    else if( aeta < 2.3   ) EA = 0.0184;
    else if( aeta < 2.4   ) EA = 0.0284;
    else if( aeta > 2.4   ) EA = 0.0591;
  }
  else if( version == 4 ){
    // Fall17 V2: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
    if(      aeta < 1.0   ) EA = 0.0668;
    else if( aeta < 1.479 ) EA = 0.1054;
    else if( aeta < 2.0   ) EA = 0.0786;
    else if( aeta < 2.2   ) EA = 0.0233;
    else if( aeta < 2.3   ) EA = 0.0078;
    else if( aeta < 2.4   ) EA = 0.0028;
    else if( aeta > 2.4   ) EA = 0.0137;
  }

  return EA;
}

float photon_EMEA03( int photonIdx, int version )
{
  float aeta = fabs(cms3.photons_p4().at(photonIdx).eta());
  float EA = -999;

  if(      version <= 1 ){
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2Archive#Selection_implementation_details
    if(      aeta < 1.0   ) EA = 0.0725;
    else if( aeta < 1.479 ) EA = 0.0604;
    else if( aeta < 2.0   ) EA = 0.0320;
    else if( aeta < 2.2   ) EA = 0.0512;
    else if( aeta < 2.3   ) EA = 0.0766;
    else if( aeta < 2.4   ) EA = 0.0949;
    else if( aeta > 2.4   ) EA = 0.1160;
  }
  else if( version == 2 ){
    // Spring16: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
    if(      aeta < 1.0   ) EA = 0.1210;
    else if( aeta < 1.479 ) EA = 0.1107;
    else if( aeta < 2.0   ) EA = 0.0699;
    else if( aeta < 2.2   ) EA = 0.1056;
    else if( aeta < 2.3   ) EA = 0.1457;
    else if( aeta < 2.4   ) EA = 0.1719;
    else if( aeta > 2.4   ) EA = 0.1998;
  }
  else if( version == 4 ){
    // Fall17 V2: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
    if(      aeta < 1.0   ) EA = 0.1113;
    else if( aeta < 1.479 ) EA = 0.0953;
    else if( aeta < 2.0   ) EA = 0.0619;
    else if( aeta < 2.2   ) EA = 0.0837;
    else if( aeta < 2.3   ) EA = 0.1070;
    else if( aeta < 2.4   ) EA = 0.1212;
    else if( aeta > 2.4   ) EA = 0.1466;
  }

  return EA;
}

float photonCHIso03EA( int photonIdx, int eaversion )
{
  float chiso = photons_recoChargedHadronIso().at(photonIdx);
  float ea    = photon_CHEA03(photonIdx, eaversion);
  float CHIso = std::max(float(0.0), chiso - evt_fixgridfastjet_all_rho() * ea);
  return CHIso;
}

float photonNHIso03EA( int photonIdx, int eaversion )
{
  float nhiso = photons_recoNeutralHadronIso().at(photonIdx);
  float ea    = photon_NHEA03(photonIdx, eaversion);
  float NHIso = std::max(float(0.0), nhiso - evt_fixgridfastjet_all_rho() * ea);
  return NHIso;
}

float photonEMIso03EA( int photonIdx, int eaversion )
{
  float emiso = photons_recoPhotonIso().at(photonIdx);
  float ea    = photon_EMEA03(photonIdx, eaversion);
  float EMIso = std::max(float(0.0), emiso - evt_fixgridfastjet_all_rho() * ea);
  return EMIso;
}

float photonEcalpfClusterIso03EA( int photonIdx )
{
  float emiso = cms3.photons_ecalPFClusterIso().at(photonIdx);
  float ea    = photon_EMEA03(photonIdx);
  float EMIso = std::max(float(0.0), emiso - evt_fixgridfastjet_all_rho() * ea);
  return EMIso;
}

float photonHcalpfClusterIso03EA( int photonIdx )
{
  float hciso = cms3.photons_hcalPFClusterIso().at(photonIdx);
  float ea    = photon_NHEA03(photonIdx);
  float HCIso = std::max(float(0.0), hciso - evt_fixgridfastjet_all_rho() * ea);
  return HCIso;
}

