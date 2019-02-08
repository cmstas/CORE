#include "JetSelections.h"
#include "ElectronSelections.h"
#include "MuonSelections.h" 
#include "Math/VectorUtil.h"
#include "Tools/utils.h"

using namespace tas;

// NOTE This is the JetIDLepVeto version for 2018 (See later JetID version that was proposed as an intermediate solution)
// More info: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
bool isTightPFJet_2018_v2(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());
  float pfjet_mf_ = pfjets_muonE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());

  if (pfjet_eta <= 2.6) {
      if (pfjet_nhf_ >= 0.90) return false;
      if (pfjet_nef_ >= 0.90) return false;
      if (pfjets_npfcands()[pfJetIdx] <= 1) return false;
      if (pfjet_mf_ >= 0.8) return false;
      if (pfjet_chf_ <= 1e-6) return false;
      if (pfjets_chargedMultiplicity()[pfJetIdx] == 0) return false;
      if (pfjet_cef_ >= 0.8) return false;
  }
  else if (pfjet_eta > 2.6 && pfjet_eta <= 2.7){
      if (pfjet_nhf_ >= 0.90) return false;
      if (pfjet_nef_ >= 0.99) return false;
      if (pfjet_mf_ >= 0.8) return false;
      if (pfjets_chargedMultiplicity()[pfJetIdx] == 0) return false;
      if (pfjet_cef_ >= 0.8) return false;
  }
  else if (pfjet_eta > 2.7 && pfjet_eta <= 3.0){
      if (pfjet_nef_ <= 0.02 || pfjet_nef_ >= 0.99) return false;
      if (pfjets_neutralMultiplicity()[pfJetIdx] <=2) return false;
  }
  else if (pfjet_eta > 3.0){
      if (pfjet_nhf_ <= 0.02) return false;
      if (pfjet_nef_ >= 0.90) return false;
      if (pfjets_neutralMultiplicity()[pfJetIdx] <= 10) return false;
  }

  return true;
}

// NOTE Preliminary! This is the JetID version (not JetIDLepVeto)
// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
bool isTightPFJet_2018_v1(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 2.6) {
      if (pfjet_nhf_ >= 0.90) return false;
      if (pfjet_nef_ >= 0.90) return false;
      if (pfjets_npfcands()[pfJetIdx] <= 1) return false;
      if (pfjet_chf_ <= 1e-6) return false;
      if (pfjets_chargedMultiplicity()[pfJetIdx] == 0) return false;
  }
  else if (pfjet_eta > 2.6 && pfjet_eta <= 2.7){
      if (pfjet_nhf_ >= 0.90) return false;
      if (pfjet_nef_ >= 0.99) return false;
      if (pfjets_chargedMultiplicity()[pfJetIdx] == 0) return false;
  }
  else if (pfjet_eta > 2.7 && pfjet_eta <= 3.0){
      if (pfjet_nef_ <= 0.02 || pfjet_nef_ >= 0.99) return false;
      if (pfjets_neutralMultiplicity()[pfJetIdx] <=2) return false;
  }
  else if (pfjet_eta > 3.0){
      if (pfjet_nhf_ <= 0.02) return false;
      if (pfjet_nef_ >= 0.90) return false;
      if (pfjets_neutralMultiplicity()[pfJetIdx] <= 10) return false;
  }

  return true;
}

// recommended jet ID for 94x 2017 analyses
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2017
bool isTightPFJet_2017_v1(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());
  

  if (pfjet_eta <= 2.4){
      if (pfjet_chf_ == 0.00) return false;
      if (pfjets_chargedMultiplicity()[pfJetIdx] == 0) return false;
  }
  if (pfjet_eta <= 2.7){
      if (pfjet_nhf_ >= 0.90) return false;
      if (pfjet_nef_ >= 0.90) return false;
      if (pfjets_npfcands()[pfJetIdx] <= 1) return false;
  }
  if (pfjet_eta > 2.7 && pfjet_eta <= 3.0){
      if (pfjet_nef_ <= 0.02 || pfjet_nef_ >= 0.99) return false;
      if (pfjets_neutralMultiplicity()[pfJetIdx] <=2) return false;
  }
  if (pfjet_eta > 3.0){
      if (pfjet_nhf_ <= 0.02) return false;
      if (pfjet_nef_ >= 0.90) return false;
      if (pfjets_neutralMultiplicity()[pfJetIdx] <= 10) return false;
  }

  return true;
}

// Jet ID as of Sept. 11th based on this twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetID?rev=95
bool isLoosePFJet_Summer16_v1(unsigned int pfJetIdx, bool use_puppi){

  float pfjet_chf_  = !use_puppi ? pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_chargedHadronE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = !use_puppi ? pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_neutralHadronE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  float pfjet_cef_  = !use_puppi ? pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_chargedEmE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  float pfjet_nef_  = !use_puppi ? pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_neutralEmE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  int   pfjet_cm_  = !use_puppi ? pfjets_chargedMultiplicity()[pfJetIdx] : pfjets_puppi_chargedMultiplicity()[pfJetIdx];
  int   pfjet_nm_  = !use_puppi ? pfjets_neutralMultiplicity()[pfJetIdx] : pfjets_puppi_neutralMultiplicity()[pfJetIdx];
  float pfjet_eta  = !use_puppi ? fabs(pfjets_p4()[pfJetIdx].eta()) : fabs(pfjets_puppi_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 2.7){
	if (pfjet_nhf_ >= 0.99       ) return false;
	if (pfjet_nef_ >= 0.99       ) return false;
	if (pfjet_cm_ + pfjet_nm_ < 2) return false;
	// if (pfjet_muf_ >= 0.8        ) return false; removed again

	if (pfjet_eta < 2.4){
	  if (!(pfjet_cm_  >   0.  ) ) return false;
	  if (!(pfjet_chf_ >   0.  ) ) return false;
	  if (!(pfjet_cef_ <   0.99) ) return false;
	}
  }else if( pfjet_eta > 2.7 && pfjet_eta <= 3.0 ){
	if (!(pfjet_nef_ < 0.9 ) ) return false;
	if (!(pfjet_nm_  > 2   ) ) return false;
  }else if( pfjet_eta > 3.0 ){
	if (!(pfjet_nef_ < 0.9 ) ) return false;
	if (!(pfjet_nm_  > 10  ) ) return false;
  }

  return true;
}

bool isTightPFJet_Summer16_v1(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 2.7){
	if (pfjet_nef_ >= 0.90) return false;
	if (pfjet_nhf_ >= 0.90) return false;
  }

  if (!isLoosePFJet_Summer16_v1(pfJetIdx)) return false;

  return true;
}

bool isTightPFJetLepVeto_Summer16_v1(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_muf_  = pfjets_muonE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 2.7){
	if (pfjet_nef_ >= 0.90) return false;
	if (pfjet_nhf_ >= 0.90) return false;
	if (pfjet_eta <= 2.4){
	  if (pfjet_cef_ >= 0.90) return false;
	}
	if (pfjet_muf_ >= 0.8        ) return false;
  }

  if (!isLoosePFJet_Summer16_v1(pfJetIdx)) return false;

  return true;
}

bool isLoosePFJet(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  int   pfjet_cm_  = pfjets_chargedMultiplicity()[pfJetIdx];
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjets_npfcands()[pfJetIdx] < 2) return false;
  if (pfjet_nef_ >= 0.99) return false;
  if (pfjet_nhf_ >= 0.99) return false;

  if (pfjet_eta < 2.4){
    if (pfjet_cm_ < 1) return false;
    if (pfjet_chf_ < 1e-6) return false;
    if (pfjet_cef_ >= 0.99) return false;
  }

  return true;
}

bool isMediumPFJet(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());

  if (pfjet_nef_ >= 0.95) return false;
  if (pfjet_nhf_ >= 0.95) return false;

  if (!isLoosePFJet(pfJetIdx)) return false;

  return true;
}

bool isTightPFJet(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());

  if (pfjet_nef_ >= 0.90) return false;
  if (pfjet_nhf_ >= 0.90) return false;

  if (!isLoosePFJet(pfJetIdx)) return false;

  return true;
}

bool isLoosePFJetV2(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_muf_  = pfjets_muonE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  int   pfjet_cm_  = pfjets_chargedMultiplicity()[pfJetIdx];
  int   pfjet_nm_  = pfjets_neutralMultiplicity()[pfJetIdx];
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_cm_ + pfjet_nm_ < 2) return false;
  if (pfjet_nef_ >= 0.99) return false;
  if (pfjet_nhf_ >= 0.99) return false;
  if (pfjet_muf_ >= 0.8) return false;

  if (pfjet_eta < 2.4){
    if (pfjet_cm_ < 1) return false;
    if (!(pfjet_chf_ > 0.)) return false;
    if (pfjet_cef_ >= 0.99) return false;
  }

  return true;
}

bool isTightPFJetV2(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_nef_ >= 0.90) return false;
  if (pfjet_nhf_ >= 0.90) return false;
  if (pfjet_eta < 2.4){
    if (pfjet_cef_ >= 0.90) return false;
  }


  if (!isLoosePFJetV2(pfJetIdx)) return false;

  return true;
}

bool isLoosePFJet_50nsV1(unsigned int pfJetIdx, bool use_puppi){

  float pfjet_chf_  = !use_puppi ? pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_chargedHadronE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = !use_puppi ? pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_neutralHadronE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  float pfjet_cef_  = !use_puppi ? pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_chargedEmE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  float pfjet_nef_  = !use_puppi ? pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy()) : pfjets_puppi_neutralEmE()[pfJetIdx] / (pfjets_puppi_undoJEC().at(pfJetIdx)*pfjets_puppi_p4()[pfJetIdx].energy());
  int   pfjet_cm_  = !use_puppi ? pfjets_chargedMultiplicity()[pfJetIdx] : pfjets_puppi_chargedMultiplicity()[pfJetIdx];
  int   pfjet_nm_  = !use_puppi ? pfjets_neutralMultiplicity()[pfJetIdx] : pfjets_puppi_neutralMultiplicity()[pfJetIdx];
  float pfjet_eta  = !use_puppi ? fabs(pfjets_p4()[pfJetIdx].eta()) : fabs(pfjets_puppi_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 3.0){
	if (pfjet_nhf_ >= 0.99       ) return false;
	if (pfjet_nef_ >= 0.99       ) return false;
	if (pfjet_cm_ + pfjet_nm_ < 2) return false;
	// if (pfjet_muf_ >= 0.8        ) return false; removed again

	if (pfjet_eta < 2.4){
	  if (!(pfjet_cm_  >   0.  ) ) return false;
	  if (!(pfjet_chf_ >   0.  ) ) return false;
	  if (!(pfjet_cef_ <   0.99) ) return false;
	}
  }else if( pfjet_eta > 3.0 ){
	if (!(pfjet_nef_ < 0.9 ) ) return false;
	if (!(pfjet_nm_  > 10  ) ) return false;
  }

  return true;
}

bool isTightPFJet_50nsV1(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 3.0){
	if (pfjet_nef_ >= 0.90) return false;
	if (pfjet_nhf_ >= 0.90) return false;
  }

  if (!isLoosePFJet_50nsV1(pfJetIdx)) return false;

  return true;
}

bool isTightPFJetLepVeto_50nsV1(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_muf_  = pfjets_muonE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 3.0){
	if (pfjet_nef_ >= 0.90) return false;
	if (pfjet_nhf_ >= 0.90) return false;
	if (pfjet_eta <= 2.4){
	  if (pfjet_cef_ >= 0.90) return false;
	}
	if (pfjet_muf_ >= 0.8        ) return false;
  }

  if (!isLoosePFJet_50nsV1(pfJetIdx)) return false;

  return true;
}

bool isMonoPFJet_MT2(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta > 3.0) return false;
  if (pfjet_chf_ <= 0.05) return false;
  if (pfjet_nhf_ >= 0.70) return false;
  if (pfjet_nef_ >= 0.80) return false;

  if (!isTightPFJet_50nsV1(pfJetIdx)) return false;

  return true;
}

bool isMonoPFJet_Monojet(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta > 3.0) return false;
  if (pfjet_chf_ <= 0.20) return false;
  if (pfjet_nhf_ >= 0.70) return false;
  if (pfjet_nef_ >= 0.70) return false;

  if (!isMonoPFJet_MT2(pfJetIdx)) return false;

  return true;
}

// this is very old, don't use
bool loosePileupJetId(unsigned int pfJetIdx){

  float eta = fabs(pfjets_p4().at(pfJetIdx).eta());
  float value = pfjets_pileupJetId().at(pfJetIdx);

  if( (eta >= 0   ) && (eta <= 2.5 ) && (value > -0.63) ) return true;
  if( (eta > 2.5  ) && (eta <= 2.75) && (value > -0.60) ) return true;
  if( (eta > 2.75 ) && (eta <= 3.0 ) && (value > -0.55) ) return true;
  if( (eta > 3.0  ) && (eta <= 5.2 ) && (value > -0.45) ) return true;
    
  return false;
}

// hard-code 81X MVA working points (used for 2016,17,18 so far)
// if they ever update for newer releases, have to make it configurable
bool pileupJetId(unsigned int pfJetIdx, id_level_t id_level){

  float pt = pfjets_p4().at(pfJetIdx).pt();
  float eta = fabs(pfjets_p4().at(pfJetIdx).eta());
  float value = pfjets_pileupJetId().at(pfJetIdx);

  if(id_level == PUID_loose){
      if(0 <= pt && pt <= 30) {
          if(0.00<=eta && eta<=2.50 && value>-0.97) return true;
          if(2.50<=eta && eta<=2.75 && value>-0.68) return true;
          if(2.75<=eta && eta<=3.00 && value>-0.53) return true;
          if(3.00<=eta &&              value>-0.47) return true;
      }
      if(30 <= pt){
          if(0.00<=eta && eta<=2.50 && value>-0.89) return true;
          if(2.50<=eta && eta<=2.75 && value>-0.52) return true;
          if(2.75<=eta && eta<=3.00 && value>-0.38) return true;
          if(3.00<=eta              && value>-0.30) return true;
      }
  }
  if(id_level == PUID_medium){
      if(0 <= pt && pt <= 30) {
          if(0.00<=eta && eta<=2.50 && value> 0.18) return true;
          if(2.50<=eta && eta<=2.75 && value>-0.55) return true;
          if(2.75<=eta && eta<=3.00 && value>-0.42) return true;
          if(3.00<=eta &&              value>-0.36) return true;
      }
      if(30 <= pt){
          if(0.00<=eta && eta<=2.50 && value> 0.61) return true;
          if(2.50<=eta && eta<=2.75 && value>-0.35) return true;
          if(2.75<=eta && eta<=3.00 && value>-0.23) return true;
          if(3.00<=eta              && value>-0.17) return true;
      }
  }
  if(id_level == PUID_tight){
      if(0 <= pt && pt <= 30) {
          if(0.00<=eta && eta<=2.50 && value> 0.69) return true;
          if(2.50<=eta && eta<=2.75 && value>-0.35) return true;
          if(2.75<=eta && eta<=3.00 && value>-0.26) return true;
          if(3.00<=eta &&              value>-0.21) return true;
      }
      if(30 <= pt){
          if(0.00<=eta && eta<=2.50 && value> 0.86) return true;
          if(2.50<=eta && eta<=2.75 && value>-0.10) return true;
          if(2.75<=eta && eta<=3.00 && value>-0.05) return true;
          if(3.00<=eta              && value>-0.01) return true;
      }
  }
    
  return false;
}

bool JetIsElectron(LorentzVector pfJet, id_level_t id_level, float ptcut, float deltaR){
  bool jetIsLep = false;
  for (unsigned int eidx = 0; eidx < tas::els_p4().size(); eidx++){
    LorentzVector electron = tas::els_p4().at(eidx);
    if (electron.pt() < ptcut) continue;
    if (!electronID(eidx,id_level)) continue;
    if (ROOT::Math::VectorUtil::DeltaR(pfJet, electron) > deltaR) continue;
    jetIsLep = true;
  }
  return jetIsLep;
}

bool JetIsMuon(LorentzVector pfJet, id_level_t id_level, float ptcut, float deltaR){
  bool jetIsLep = false;
  for (unsigned int muidx = 0; muidx < tas::mus_p4().size(); muidx++){
    LorentzVector muon = tas::mus_p4().at(muidx);
    if (muon.pt() < ptcut) continue;
    if (!muonID(muidx,id_level)) continue;
    if (ROOT::Math::VectorUtil::DeltaR(pfJet, muon) > deltaR) continue;
    jetIsLep = true;
  }
  return jetIsLep;
}

bool loosePileupJetId_v2(unsigned int pfJetIdx, bool use_puppi){


  float eta = use_puppi ? fabs(pfjets_puppi_p4().at(pfJetIdx).eta()) : fabs(pfjets_p4().at(pfJetIdx).eta());
  float value = use_puppi ? pfjets_puppi_pileupJetId().at(pfJetIdx) : pfjets_pileupJetId().at(pfJetIdx);
  float pt = use_puppi ? pfjets_puppi_p4().at(pfJetIdx).pt() : pfjets_p4().at(pfJetIdx).pt();

  if( (eta >= 0   ) && (eta <= 2.5 ) && (pt > 30) && (value > -0.63) ) return true;
  if( (eta > 2.5  ) && (eta <= 2.75) && (pt > 30) && (value > -0.60) ) return true;
  if( (eta > 2.75 ) && (eta <= 3.0 ) && (pt > 30) && (value > -0.55) ) return true;
  if( (eta > 3.0  ) && (eta <= 5.2 ) && (pt > 30) && (value > -0.45) ) return true;

  if( (eta >= 0   ) && (eta <= 2.5 ) && (pt > 20) && (value > -0.63) ) return true;
  if( (eta > 2.5  ) && (eta <= 2.75) && (pt > 20) && (value > -0.60) ) return true;
  if( (eta > 2.75 ) && (eta <= 3.0 ) && (pt > 20) && (value > -0.55) ) return true;
  if( (eta > 3.0  ) && (eta <= 5.2 ) && (pt > 20) && (value > -0.45) ) return true;

  if( (eta >= 0   ) && (eta <= 2.5 ) && (pt > 10) && (value > -0.95) ) return true;
  if( (eta > 2.5  ) && (eta <= 2.75) && (pt > 10) && (value > -0.96) ) return true;
  if( (eta > 2.75 ) && (eta <= 3.0 ) && (pt > 10) && (value > -0.94) ) return true;
  if( (eta > 3.0  ) && (eta <= 5.2 ) && (pt > 10) && (value > -0.95) ) return true;

  if( (eta >= 0   ) && (eta <= 2.5 ) && (pt > 0) && (value > -0.95) ) return true;
  if( (eta > 2.5  ) && (eta <= 2.75) && (pt > 0) && (value > -0.96) ) return true;
  if( (eta > 2.75 ) && (eta <= 3.0 ) && (pt > 0) && (value > -0.94) ) return true;
  if( (eta > 3.0  ) && (eta <= 5.2 ) && (pt > 0) && (value > -0.95) ) return true;
    
  return false;
}

// returns true if the PFJet is BAD -- typically also require a pt > 20 cut in the analysis
bool isBadFastsimJet(unsigned int pfJetIdx){

  float pfjet_chf  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());
  float pfjet_mcdr = pfjets_mcdr()[pfJetIdx];

  if (pfjet_eta < 2.5 && pfjet_chf < 0.1 && fabs(pfjet_mcdr) > 0.3) return true;
  return false;
}

float getPrefireInefficiency_singlejet_2017(float pt, float eta) {
  if (pt >= 30 && pt < 35 && eta >= -3.100 && eta < -3.000) return 0.02500; // +- 88.5%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.00826; // +- 28.1%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.00538; // +- 26.5%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.01634; // +- 37.0%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.00787; // +- 20.1%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.00724; // +- 13.1%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.00669; // +- 35.5%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.01070; // +- 51.5%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.00775; // +- 34.4%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.00746; // +- 18.1%
  if (pt >= 35 && pt < 40 && eta >= 3.000 && eta < 3.100) return 0.03571; // +- 91.3%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.00620; // +- 11.0%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.00921; // +- 20.3%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.01050; // +- 32.3%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.00800; // +- 35.3%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.00515; // +- 12.2%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.01600; // +- 26.1%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.00838; // +- 16.4%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.00601; // +- 8.2%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.00666; // +- 15.8%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.01705; // +- 18.2%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.01175; // +- 25.2%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.01554; // +- 27.3%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.04408; // +- 15.6%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.01898; // +- 15.8%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.00971; // +- 15.9%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.04403; // +- 17.4%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.05013; // +- 14.8%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.05583; // +- 15.9%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.01521; // +- 28.8%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.01890; // +- 26.9%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.10577; // +- 11.4%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.04289; // +- 15.6%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.04899; // +- 15.9%
  if (pt >= 100 && pt < 150 && eta >= -3.000 && eta < -2.750) return 0.22151; // +- 8.6%
  if (pt >= 100 && pt < 150 && eta >= -2.750 && eta < -2.500) return 0.15058; // +- 8.4%
  if (pt >= 100 && pt < 150 && eta >= -2.500 && eta < -2.250) return 0.10558; // +- 11.1%
  if (pt >= 100 && pt < 150 && eta >= -2.250 && eta < -2.000) return 0.02595; // +- 21.6%
  if (pt >= 100 && pt < 150 && eta >= 2.000 && eta < 2.250) return 0.04027; // +- 18.8%
  if (pt >= 100 && pt < 150 && eta >= 2.250 && eta < 2.500) return 0.18996; // +- 7.5%
  if (pt >= 100 && pt < 150 && eta >= 2.500 && eta < 2.750) return 0.13980; // +- 8.7%
  if (pt >= 100 && pt < 150 && eta >= 2.750 && eta < 3.000) return 0.25420; // +- 7.8%
  if (pt >= 150 && pt < 200 && eta >= -3.100 && eta < -3.000) return 0.25000; // +- 86.0%
  if (pt >= 150 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.45076; // +- 6.8%
  if (pt >= 150 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.30090; // +- 6.4%
  if (pt >= 150 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.14911; // +- 9.5%
  if (pt >= 150 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.02819; // +- 21.7%
  if (pt >= 150 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.04944; // +- 16.6%
  if (pt >= 150 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.22300; // +- 7.7%
  if (pt >= 150 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.32596; // +- 6.4%
  if (pt >= 150 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.45914; // +- 6.8%
  if (pt >= 200 && pt < 300 && eta >= -3.100 && eta < -3.000) return 0.66667; // +- 41.0%
  if (pt >= 200 && pt < 300 && eta >= -3.000 && eta < -2.750) return 0.59353; // +- 5.0%
  if (pt >= 200 && pt < 300 && eta >= -2.750 && eta < -2.500) return 0.44207; // +- 4.4%
  if (pt >= 200 && pt < 300 && eta >= -2.500 && eta < -2.250) return 0.18806; // +- 6.9%
  if (pt >= 200 && pt < 300 && eta >= -2.250 && eta < -2.000) return 0.04554; // +- 13.3%
  if (pt >= 200 && pt < 300 && eta >= 2.000 && eta < 2.250) return 0.08665; // +- 10.1%
  if (pt >= 200 && pt < 300 && eta >= 2.250 && eta < 2.500) return 0.30733; // +- 5.1%
  if (pt >= 200 && pt < 300 && eta >= 2.500 && eta < 2.750) return 0.50664; // +- 4.0%
  if (pt >= 200 && pt < 300 && eta >= 2.750 && eta < 3.000) return 0.66087; // +- 4.7%
  if (pt >= 200 && pt < 300 && eta >= 3.000 && eta < 3.100) return 0.33333; // +- 81.3%
  if (pt >= 300  && eta >= -3.100 && eta < -3.000) return 0.33333; // +- 81.3%
  if (pt >= 300  && eta >= -3.000 && eta < -2.750) return 0.81373; // +- 4.8%
  if (pt >= 300  && eta >= -2.750 && eta < -2.500) return 0.61630; // +- 3.5%
  if (pt >= 300  && eta >= -2.500 && eta < -2.250) return 0.27482; // +- 4.8%
  if (pt >= 300  && eta >= -2.250 && eta < -2.000) return 0.06092; // +- 8.6%
  if (pt >= 300  && eta >= 2.000 && eta < 2.250) return 0.13857; // +- 5.6%
  if (pt >= 300  && eta >= 2.250 && eta < 2.500) return 0.43550; // +- 3.3%
  if (pt >= 300  && eta >= 2.500 && eta < 2.750) return 0.74145; // +- 2.7%
  if (pt >= 300  && eta >= 2.750 && eta < 3.000) return 0.85417; // +- 4.3%
  return 0.;
}

float getPrefireInefficiencyError_singlejet_2017(float pt, float eta) {
  if (pt >= 30 && pt < 35 && eta >= -3.100 && eta < -3.000) return 0.02214; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -3.000 && eta < -2.750) return 0.00147; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.00232; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.500 && eta < -2.250) return 0.00486; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.00142; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.000 && eta < 2.250) return 0.00496; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.00605; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.00158; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.750 && eta < 3.000) return 0.00106; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 3.000 && eta < 3.100) return 0.02816; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -3.100 && eta < -3.000) return 0.02531; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.00095; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.750 && eta < -2.500) return 0.00282; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.00237; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.00550; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.000 && eta < 2.250) return 0.00618; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.00266; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.00135; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.750 && eta < 3.000) return 0.00098; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 3.000 && eta < 3.100) return 0.03261; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.100 && eta < -3.000) return 0.02298; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.00068; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.00187; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.00339; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.00282; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.00063; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.00417; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.00137; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.00049; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 3.000 && eta < 3.100) return 0.02247; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -3.100 && eta < -3.000) return 0.02197; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.00105; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.00311; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.00296; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.00316; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.00424; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.00688; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.00300; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.00155; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 3.000 && eta < 3.100) return 0.01709; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -3.100 && eta < -3.000) return 0.03388; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.00768; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.00744; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.00890; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.00438; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.00509; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.01205; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.00668; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.00778; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 3.000 && eta < 3.100) return 0.02898; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -3.100 && eta < -3.000) return 0.05399; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -3.000 && eta < -2.750) return 0.01910; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -2.750 && eta < -2.500) return 0.01265; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -2.500 && eta < -2.250) return 0.01168; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -2.250 && eta < -2.000) return 0.00560; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.000 && eta < 2.250) return 0.00756; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.250 && eta < 2.500) return 0.01430; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.500 && eta < 2.750) return 0.01212; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.750 && eta < 3.000) return 0.01983; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 3.000 && eta < 3.100) return 0.05399; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -3.100 && eta < -3.000) return 0.21504; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.03059; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.01938; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.01414; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.00613; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.00823; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.01723; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.02094; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.03106; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 3.000 && eta < 3.100) return 0.35355; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -3.100 && eta < -3.000) return 0.27317; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -3.000 && eta < -2.750) return 0.02951; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.750 && eta < -2.500) return 0.01937; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.500 && eta < -2.250) return 0.01298; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.250 && eta < -2.000) return 0.00607; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.000 && eta < 2.250) return 0.00874; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.250 && eta < 2.500) return 0.01567; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.500 && eta < 2.750) return 0.02038; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.750 && eta < 3.000) return 0.03133; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 3.000 && eta < 3.100) return 0.27113; // +- 0.0%
  if (pt >= 300  && eta >= -3.100 && eta < -3.000) return 0.27113; // +- 0.0%
  if (pt >= 300  && eta >= -3.000 && eta < -2.750) return 0.03894; // +- 0.0%
  if (pt >= 300  && eta >= -2.750 && eta < -2.500) return 0.02173; // +- 0.0%
  if (pt >= 300  && eta >= -2.500 && eta < -2.250) return 0.01322; // +- 0.0%
  if (pt >= 300  && eta >= -2.250 && eta < -2.000) return 0.00522; // +- 0.0%
  if (pt >= 300  && eta >= 2.000 && eta < 2.250) return 0.00781; // +- 0.0%
  if (pt >= 300  && eta >= 2.250 && eta < 2.500) return 0.01457; // +- 0.0%
  if (pt >= 300  && eta >= 2.500 && eta < 2.750) return 0.02036; // +- 0.0%
  if (pt >= 300  && eta >= 2.750 && eta < 3.000) return 0.03653; // +- 0.0%
  if (pt >= 300  && eta >= 3.000 && eta < 3.100) return 0.35355; // +- 0.0%
  return 0.;
}

float getPrefireInefficiency_singlephoton_2017(float pt, float eta) {
  if (pt >= 20 && pt < 25 && eta >= -3.000 && eta < -2.750) return 0.01333; // +- 45.5%
  if (pt >= 20 && pt < 25 && eta >= -2.750 && eta < -2.500) return 0.01770; // +- 42.1%
  if (pt >= 20 && pt < 25 && eta >= -2.500 && eta < -2.250) return 0.03982; // +- 30.6%
  if (pt >= 20 && pt < 25 && eta >= 2.000 && eta < 2.250) return 0.01020; // +- 41.1%
  if (pt >= 20 && pt < 25 && eta >= 2.250 && eta < 2.500) return 0.03309; // +- 30.3%
  if (pt >= 20 && pt < 25 && eta >= 2.500 && eta < 2.750) return 0.01521; // +- 40.8%
  if (pt >= 20 && pt < 25 && eta >= 2.750 && eta < 3.000) return 0.01802; // +- 42.2%
  if (pt >= 25 && pt < 30 && eta >= -3.000 && eta < -2.750) return 0.01863; // +- 49.0%
  if (pt >= 25 && pt < 30 && eta >= -2.750 && eta < -2.500) return 0.03883; // +- 32.4%
  if (pt >= 25 && pt < 30 && eta >= -2.500 && eta < -2.250) return 0.05936; // +- 25.8%
  if (pt >= 25 && pt < 30 && eta >= -2.250 && eta < -2.000) return 0.03030; // +- 31.9%
  if (pt >= 25 && pt < 30 && eta >= 2.000 && eta < 2.250) return 0.02479; // +- 36.1%
  if (pt >= 25 && pt < 30 && eta >= 2.250 && eta < 2.500) return 0.13216; // +- 16.7%
  if (pt >= 25 && pt < 30 && eta >= 2.500 && eta < 2.750) return 0.05000; // +- 29.3%
  if (pt >= 25 && pt < 30 && eta >= 2.750 && eta < 3.000) return 0.03185; // +- 40.5%
  if (pt >= 30 && pt < 35 && eta >= -3.000 && eta < -2.750) return 0.03226; // +- 45.3%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.04369; // +- 30.8%
  if (pt >= 30 && pt < 35 && eta >= -2.500 && eta < -2.250) return 0.11278; // +- 16.9%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.02980; // +- 30.0%
  if (pt >= 30 && pt < 35 && eta >= 2.000 && eta < 2.250) return 0.03550; // +- 26.3%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.22179; // +- 11.6%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.07853; // +- 24.0%
  if (pt >= 30 && pt < 35 && eta >= 2.750 && eta < 3.000) return 0.04878; // +- 37.8%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.15534; // +- 22.7%
  if (pt >= 35 && pt < 40 && eta >= -2.750 && eta < -2.500) return 0.05714; // +- 26.8%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.12422; // +- 14.5%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.02500; // +- 26.7%
  if (pt >= 35 && pt < 40 && eta >= 2.000 && eta < 2.250) return 0.05825; // +- 19.0%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.24918; // +- 9.9%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.07643; // +- 26.9%
  if (pt >= 35 && pt < 40 && eta >= 2.750 && eta < 3.000) return 0.12821; // +- 23.7%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.25568; // +- 12.8%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.14121; // +- 13.0%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.14157; // +- 11.5%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.03494; // +- 19.9%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.07120; // +- 14.0%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.27157; // +- 8.2%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.15599; // +- 12.1%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.24309; // +- 13.0%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.44140; // +- 5.6%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.27212; // +- 7.6%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.13988; // +- 13.3%
  if (pt >= 50 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.01806; // +- 29.9%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.07767; // +- 16.5%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.23901; // +- 9.3%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.31876; // +- 6.2%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.41755; // +- 6.1%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.57500; // +- 5.6%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.43728; // +- 6.8%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.21795; // +- 12.3%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.03289; // +- 28.7%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.09764; // +- 17.2%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.32273; // +- 9.7%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.52727; // +- 5.2%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.60079; // +- 5.1%
  if (pt >= 100 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.78967; // +- 3.2%
  if (pt >= 100 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.55118; // +- 4.6%
  if (pt >= 100 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.21221; // +- 10.3%
  if (pt >= 100 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.03542; // +- 22.1%
  if (pt >= 100 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.10484; // +- 12.8%
  if (pt >= 100 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.38338; // +- 6.6%
  if (pt >= 100 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.69762; // +- 3.1%
  if (pt >= 100 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.74460; // +- 3.5%
  if (pt >= 200  && eta >= -3.000 && eta < -2.750) return 0.86792; // +- 3.1%
  if (pt >= 200  && eta >= -2.750 && eta < -2.500) return 0.70561; // +- 4.4%
  if (pt >= 200  && eta >= -2.500 && eta < -2.250) return 0.31298; // +- 9.1%
  if (pt >= 200  && eta >= -2.250 && eta < -2.000) return 0.04922; // +- 21.3%
  if (pt >= 200  && eta >= 2.000 && eta < 2.250) return 0.13449; // +- 11.6%
  if (pt >= 200  && eta >= 2.250 && eta < 2.500) return 0.50579; // +- 6.1%
  if (pt >= 200  && eta >= 2.500 && eta < 2.750) return 0.79649; // +- 3.0%
  if (pt >= 200  && eta >= 2.750 && eta < 3.000) return 0.76159; // +- 4.6%
  return 0.;
}

float getPrefireInefficiencyError_singlephoton_2017(float pt, float eta) {
  if (pt >= 20 && pt < 25 && eta >= -3.000 && eta < -2.750) return 0.00606; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= -2.750 && eta < -2.500) return 0.00745; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= -2.500 && eta < -2.250) return 0.01220; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= -2.250 && eta < -2.000) return 0.00357; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.000 && eta < 2.250) return 0.00420; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.250 && eta < 2.500) return 0.01002; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.500 && eta < 2.750) return 0.00620; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.750 && eta < 3.000) return 0.00761; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -3.000 && eta < -2.750) return 0.00914; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -2.750 && eta < -2.500) return 0.01260; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -2.500 && eta < -2.250) return 0.01532; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -2.250 && eta < -2.000) return 0.00967; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.000 && eta < 2.250) return 0.00895; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.250 && eta < 2.500) return 0.02211; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.500 && eta < 2.750) return 0.01466; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.750 && eta < 3.000) return 0.01290; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -3.000 && eta < -2.750) return 0.01462; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.01344; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.500 && eta < -2.250) return 0.01901; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.00895; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.000 && eta < 2.250) return 0.00935; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.02570; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.01889; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.750 && eta < 3.000) return 0.01845; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.03522; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.750 && eta < -2.500) return 0.01534; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.01806; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.00667; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.000 && eta < 2.250) return 0.01106; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.02460; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.02055; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.750 && eta < 3.000) return 0.03039; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.03267; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.01841; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.01628; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.00695; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.01000; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.02228; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.01890; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.03166; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.02477; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.02081; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.01864; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.00539; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.01279; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.02219; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.01980; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.02539; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.03196; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.02966; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.02676; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.00944; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.01682; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.03139; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.02750; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.03085; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.02497; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.02551; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.02185; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.00784; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.01346; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.02511; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.02144; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.02632; // +- 0.0%
  if (pt >= 200  && eta >= -3.000 && eta < -2.750) return 0.02728; // +- 0.0%
  if (pt >= 200  && eta >= -2.750 && eta < -2.500) return 0.03131; // +- 0.0%
  if (pt >= 200  && eta >= -2.500 && eta < -2.250) return 0.02852; // +- 0.0%
  if (pt >= 200  && eta >= -2.250 && eta < -2.000) return 0.01046; // +- 0.0%
  if (pt >= 200  && eta >= 2.000 && eta < 2.250) return 0.01564; // +- 0.0%
  if (pt >= 200  && eta >= 2.250 && eta < 2.500) return 0.03107; // +- 0.0%
  if (pt >= 200  && eta >= 2.500 && eta < 2.750) return 0.02406; // +- 0.0%
  if (pt >= 200  && eta >= 2.750 && eta < 3.000) return 0.03492; // +- 0.0%
  return 0.;
}

float getPrefireInefficiency_singlejet_2016(float pt, float eta) {
  if (pt >= 30 && pt < 35 && eta >= -3.000 && eta < -2.750) return 0.01332; // +- 19.1%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.01790; // +- 18.0%
  if (pt >= 30 && pt < 35 && eta >= -2.500 && eta < -2.250) return 0.00826; // +- 25.6%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.01648; // +- 33.9%
  if (pt >= 30 && pt < 35 && eta >= 2.000 && eta < 2.250) return 0.01155; // +- 33.6%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.00834; // +- 23.9%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.00720; // +- 16.7%
  if (pt >= 30 && pt < 35 && eta >= 2.750 && eta < 3.000) return 0.01397; // +- 17.4%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.01277; // +- 20.1%
  if (pt >= 35 && pt < 40 && eta >= -2.750 && eta < -2.500) return 0.01641; // +- 19.5%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.01475; // +- 25.6%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.01902; // +- 32.2%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.01825; // +- 22.6%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.01074; // +- 17.7%
  if (pt >= 35 && pt < 40 && eta >= 2.750 && eta < 3.000) return 0.00644; // +- 15.7%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.01751; // +- 16.5%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.01209; // +- 16.7%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.01601; // +- 20.0%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.01095; // +- 27.8%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.01422; // +- 25.3%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.01355; // +- 19.8%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.01615; // +- 15.3%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.00892; // +- 17.1%
  if (pt >= 50 && pt < 70 && eta >= -3.100 && eta < -3.000) return 0.02941; // +- 90.0%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.01082; // +- 20.3%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.01589; // +- 16.8%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.01967; // +- 17.9%
  if (pt >= 50 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.01094; // +- 23.2%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.01693; // +- 20.2%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.01761; // +- 17.5%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.01637; // +- 15.4%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.00857; // +- 19.4%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.02308; // +- 25.3%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.02328; // +- 19.2%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.02676; // +- 19.0%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.01273; // +- 25.9%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.01910; // +- 22.8%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.04369; // +- 15.4%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.02682; // +- 18.2%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.02448; // +- 23.6%
  if (pt >= 100 && pt < 150 && eta >= -3.000 && eta < -2.750) return 0.08028; // +- 15.7%
  if (pt >= 100 && pt < 150 && eta >= -2.750 && eta < -2.500) return 0.07178; // +- 12.1%
  if (pt >= 100 && pt < 150 && eta >= -2.500 && eta < -2.250) return 0.03694; // +- 17.0%
  if (pt >= 100 && pt < 150 && eta >= -2.250 && eta < -2.000) return 0.01223; // +- 25.5%
  if (pt >= 100 && pt < 150 && eta >= 2.000 && eta < 2.250) return 0.02796; // +- 19.5%
  if (pt >= 100 && pt < 150 && eta >= 2.250 && eta < 2.500) return 0.07925; // +- 11.3%
  if (pt >= 100 && pt < 150 && eta >= 2.500 && eta < 2.750) return 0.06612; // +- 12.4%
  if (pt >= 100 && pt < 150 && eta >= 2.750 && eta < 3.000) return 0.08597; // +- 15.1%
  if (pt >= 100 && pt < 150 && eta >= 3.000 && eta < 3.100) return 0.09091; // +- 92.9%
  if (pt >= 150 && pt < 200 && eta >= -3.100 && eta < -3.000) return 0.25000; // +- 86.0%
  if (pt >= 150 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.26587; // +- 10.4%
  if (pt >= 150 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.17434; // +- 8.7%
  if (pt >= 150 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.07034; // +- 13.7%
  if (pt >= 150 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.01597; // +- 26.1%
  if (pt >= 150 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.02203; // +- 22.5%
  if (pt >= 150 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.14630; // +- 9.2%
  if (pt >= 150 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.14311; // +- 10.1%
  if (pt >= 150 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.21918; // +- 12.6%
  if (pt >= 200 && pt < 300 && eta >= -3.000 && eta < -2.750) return 0.39676; // +- 7.8%
  if (pt >= 200 && pt < 300 && eta >= -2.750 && eta < -2.500) return 0.26440; // +- 6.4%
  if (pt >= 200 && pt < 300 && eta >= -2.500 && eta < -2.250) return 0.09872; // +- 9.3%
  if (pt >= 200 && pt < 300 && eta >= -2.250 && eta < -2.000) return 0.02222; // +- 16.5%
  if (pt >= 200 && pt < 300 && eta >= 2.000 && eta < 2.250) return 0.03586; // +- 13.8%
  if (pt >= 200 && pt < 300 && eta >= 2.250 && eta < 2.500) return 0.20344; // +- 6.2%
  if (pt >= 200 && pt < 300 && eta >= 2.500 && eta < 2.750) return 0.28917; // +- 5.9%
  if (pt >= 200 && pt < 300 && eta >= 2.750 && eta < 3.000) return 0.34884; // +- 8.5%
  if (pt >= 300  && eta >= -3.000 && eta < -2.750) return 0.53788; // +- 8.1%
  if (pt >= 300  && eta >= -2.750 && eta < -2.500) return 0.48093; // +- 3.8%
  if (pt >= 300  && eta >= -2.500 && eta < -2.250) return 0.13445; // +- 6.2%
  if (pt >= 300  && eta >= -2.250 && eta < -2.000) return 0.02536; // +- 10.8%
  if (pt >= 300  && eta >= 2.000 && eta < 2.250) return 0.05669; // +- 7.3%
  if (pt >= 300  && eta >= 2.250 && eta < 2.500) return 0.26204; // +- 4.2%
  if (pt >= 300  && eta >= 2.500 && eta < 2.750) return 0.44367; // +- 4.2%
  if (pt >= 300  && eta >= 2.750 && eta < 3.000) return 0.51754; // +- 9.0%
  return 0.;
}

float getPrefireInefficiencyError_singlejet_2016(float pt, float eta) {
  if (pt >= 30 && pt < 35 && eta >= -3.100 && eta < -3.000) return 0.03388; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -3.000 && eta < -2.750) return 0.00254; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.00322; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.500 && eta < -2.250) return 0.00212; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.00558; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.000 && eta < 2.250) return 0.00388; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.00199; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.00120; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.750 && eta < 3.000) return 0.00243; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 3.000 && eta < 3.100) return 0.02898; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -3.100 && eta < -3.000) return 0.04441; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.00256; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.750 && eta < -2.500) return 0.00321; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.00377; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.00613; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.000 && eta < 2.250) return 0.00347; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.00413; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.00190; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.750 && eta < 3.000) return 0.00101; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 3.000 && eta < 3.100) return 0.03507; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.100 && eta < -3.000) return 0.03388; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.00288; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.00201; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.00320; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.00304; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.00361; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.00268; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.00248; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.00152; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 3.000 && eta < 3.100) return 0.02597; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -3.100 && eta < -3.000) return 0.02647; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.00220; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.00267; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.00352; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.00254; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.00343; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.00309; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.00252; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.00166; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 3.000 && eta < 3.100) return 0.02469; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -3.100 && eta < -3.000) return 0.06883; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.00584; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.00446; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.00509; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.00329; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.00435; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.00672; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.00488; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.00578; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 3.000 && eta < 3.100) return 0.03771; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -3.100 && eta < -3.000) return 0.10476; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -3.000 && eta < -2.750) return 0.01264; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -2.750 && eta < -2.500) return 0.00871; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -2.500 && eta < -2.250) return 0.00628; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= -2.250 && eta < -2.000) return 0.00312; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.000 && eta < 2.250) return 0.00547; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.250 && eta < 2.500) return 0.00895; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.500 && eta < 2.750) return 0.00823; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 2.750 && eta < 3.000) return 0.01298; // +- 0.0%
  if (pt >= 100 && pt < 150 && eta >= 3.000 && eta < 3.100) return 0.08449; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -3.100 && eta < -3.000) return 0.21504; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.02766; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.01521; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.00966; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.00416; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.00496; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.01353; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.01450; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.02772; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 3.000 && eta < 3.100) return 0.13226; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -3.100 && eta < -3.000) return 0.27217; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -3.000 && eta < -2.750) return 0.03106; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.750 && eta < -2.500) return 0.01685; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.500 && eta < -2.250) return 0.00916; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.250 && eta < -2.000) return 0.00367; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.000 && eta < 2.250) return 0.00494; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.250 && eta < 2.500) return 0.01269; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.500 && eta < 2.750) return 0.01702; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.750 && eta < 3.000) return 0.02957; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 3.000 && eta < 3.100) return 0.21651; // +- 0.0%
  if (pt >= 300  && eta >= -3.000 && eta < -2.750) return 0.04343; // +- 0.0%
  if (pt >= 300  && eta >= -2.750 && eta < -2.500) return 0.01843; // +- 0.0%
  if (pt >= 300  && eta >= -2.500 && eta < -2.250) return 0.00836; // +- 0.0%
  if (pt >= 300  && eta >= -2.250 && eta < -2.000) return 0.00273; // +- 0.0%
  if (pt >= 300  && eta >= 2.000 && eta < 2.250) return 0.00415; // +- 0.0%
  if (pt >= 300  && eta >= 2.250 && eta < 2.500) return 0.01107; // +- 0.0%
  if (pt >= 300  && eta >= 2.500 && eta < 2.750) return 0.01851; // +- 0.0%
  if (pt >= 300  && eta >= 2.750 && eta < 3.000) return 0.04681; // +- 0.0%
  if (pt >= 300  && eta >= 3.000 && eta < 3.100) return 0.35355; // +- 0.0%
  return 0.;
}

float getPrefireInefficiency_singlephoton_2016(float pt, float eta) {
  if (pt >= 20 && pt < 25 && eta >= -3.000 && eta < -2.750) return 0.01154; // +- 43.3%
  if (pt >= 20 && pt < 25 && eta >= -2.750 && eta < -2.500) return 0.00990; // +- 40.5%
  if (pt >= 20 && pt < 25 && eta >= -2.500 && eta < -2.250) return 0.01597; // +- 36.9%
  if (pt >= 20 && pt < 25 && eta >= -2.250 && eta < -2.000) return 0.00940; // +- 39.4%
  if (pt >= 20 && pt < 25 && eta >= 2.000 && eta < 2.250) return 0.01117; // +- 37.1%
  if (pt >= 20 && pt < 25 && eta >= 2.250 && eta < 2.500) return 0.02007; // +- 35.1%
  if (pt >= 25 && pt < 30 && eta >= -3.000 && eta < -2.750) return 0.03279; // +- 37.1%
  if (pt >= 25 && pt < 30 && eta >= -2.750 && eta < -2.500) return 0.02013; // +- 35.1%
  if (pt >= 25 && pt < 30 && eta >= -2.500 && eta < -2.250) return 0.03427; // +- 27.5%
  if (pt >= 25 && pt < 30 && eta >= 2.250 && eta < 2.500) return 0.05112; // +- 23.2%
  if (pt >= 25 && pt < 30 && eta >= 2.500 && eta < 2.750) return 0.03887; // +- 27.7%
  if (pt >= 25 && pt < 30 && eta >= 2.750 && eta < 3.000) return 0.00541; // +- 27.4%
  if (pt >= 30 && pt < 35 && eta >= -3.000 && eta < -2.750) return 0.01325; // +- 55.6%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.03010; // +- 30.1%
  if (pt >= 30 && pt < 35 && eta >= -2.500 && eta < -2.250) return 0.04187; // +- 22.3%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.00952; // +- 34.4%
  if (pt >= 30 && pt < 35 && eta >= 2.000 && eta < 2.250) return 0.02459; // +- 25.5%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.10326; // +- 15.0%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.01954; // +- 35.0%
  if (pt >= 30 && pt < 35 && eta >= 2.750 && eta < 3.000) return 0.01183; // +- 53.6%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.02186; // +- 43.5%
  if (pt >= 35 && pt < 40 && eta >= -2.750 && eta < -2.500) return 0.04389; // +- 24.7%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.05843; // +- 18.2%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.01556; // +- 29.0%
  if (pt >= 35 && pt < 40 && eta >= 2.000 && eta < 2.250) return 0.03422; // +- 21.5%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.10174; // +- 14.5%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.04590; // +- 24.7%
  if (pt >= 35 && pt < 40 && eta >= 2.750 && eta < 3.000) return 0.03145; // +- 40.5%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.09690; // +- 18.6%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.06346; // +- 17.3%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.06545; // +- 15.5%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.01107; // +- 26.1%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.03421; // +- 19.4%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.16334; // +- 9.5%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.07342; // +- 17.3%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.06920; // +- 20.8%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.15361; // +- 12.7%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.18907; // +- 9.8%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.09375; // +- 15.5%
  if (pt >= 50 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.01174; // +- 33.8%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.04555; // +- 20.2%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.16102; // +- 12.0%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.17078; // +- 9.9%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.18750; // +- 10.7%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.35754; // +- 10.0%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.25103; // +- 11.0%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.07826; // +- 22.0%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.02574; // +- 33.6%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.06304; // +- 19.9%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.20385; // +- 12.1%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.30693; // +- 8.6%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.30303; // +- 10.7%
  if (pt >= 100 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.65500; // +- 5.1%
  if (pt >= 100 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.45330; // +- 5.8%
  if (pt >= 100 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.15691; // +- 11.8%
  if (pt >= 100 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.02326; // +- 23.5%
  if (pt >= 100 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.04847; // +- 16.9%
  if (pt >= 100 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.24800; // +- 8.9%
  if (pt >= 100 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.48969; // +- 5.2%
  if (pt >= 100 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.51485; // +- 6.8%
  if (pt >= 200  && eta >= -3.000 && eta < -2.750) return 0.73494; // +- 6.6%
  if (pt >= 200  && eta >= -2.750 && eta < -2.500) return 0.56863; // +- 7.1%
  if (pt >= 200  && eta >= -2.500 && eta < -2.250) return 0.17447; // +- 14.0%
  if (pt >= 200  && eta >= -2.250 && eta < -2.000) return 0.03879; // +- 21.6%
  if (pt >= 200  && eta >= 2.000 && eta < 2.250) return 0.06485; // +- 16.7%
  if (pt >= 200  && eta >= 2.250 && eta < 2.500) return 0.30739; // +- 9.3%
  if (pt >= 200  && eta >= 2.500 && eta < 2.750) return 0.68394; // +- 4.9%
  if (pt >= 200  && eta >= 2.750 && eta < 3.000) return 0.54667; // +- 10.5%
  return 0.;
}

float getPrefireInefficiencyError_singlephoton_2016(float pt, float eta) {
  if (pt >= 20 && pt < 25 && eta >= -3.000 && eta < -2.750) return 0.00500; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= -2.750 && eta < -2.500) return 0.00401; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= -2.500 && eta < -2.250) return 0.00589; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= -2.250 && eta < -2.000) return 0.00371; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.000 && eta < 2.250) return 0.00414; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.250 && eta < 2.500) return 0.00704; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.500 && eta < 2.750) return 0.00382; // +- 0.0%
  if (pt >= 20 && pt < 25 && eta >= 2.750 && eta < 3.000) return 0.00577; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -3.000 && eta < -2.750) return 0.01215; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -2.750 && eta < -2.500) return 0.00707; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -2.500 && eta < -2.250) return 0.00941; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= -2.250 && eta < -2.000) return 0.00369; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.000 && eta < 2.250) return 0.00336; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.250 && eta < 2.500) return 0.01186; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.500 && eta < 2.750) return 0.01075; // +- 0.0%
  if (pt >= 25 && pt < 30 && eta >= 2.750 && eta < 3.000) return 0.00148; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -3.000 && eta < -2.750) return 0.00736; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.750 && eta < -2.500) return 0.00905; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.500 && eta < -2.250) return 0.00935; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= -2.250 && eta < -2.000) return 0.00327; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.000 && eta < 2.250) return 0.00627; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.250 && eta < 2.500) return 0.01552; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.500 && eta < 2.750) return 0.00683; // +- 0.0%
  if (pt >= 30 && pt < 35 && eta >= 2.750 && eta < 3.000) return 0.00634; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -3.000 && eta < -2.750) return 0.00952; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.750 && eta < -2.500) return 0.01082; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.500 && eta < -2.250) return 0.01066; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= -2.250 && eta < -2.000) return 0.00451; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.000 && eta < 2.250) return 0.00734; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.250 && eta < 2.500) return 0.01472; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.500 && eta < 2.750) return 0.01134; // +- 0.0%
  if (pt >= 35 && pt < 40 && eta >= 2.750 && eta < 3.000) return 0.01273; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.01799; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.01097; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.01016; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.00289; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.00664; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.01555; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.01270; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.01442; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.01952; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.01850; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.01451; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.00396; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.00919; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.01929; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.01687; // +- 0.0%
  if (pt >= 50 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.02013; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.03571; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.02763; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.01718; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.00864; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.01252; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.02475; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.02637; // +- 0.0%
  if (pt >= 70 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.03251; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.03373; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.02607; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.01851; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.00546; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.00820; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.02215; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.02537; // +- 0.0%
  if (pt >= 100 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.03517; // +- 0.0%
  if (pt >= 200  && eta >= -3.000 && eta < -2.750) return 0.04873; // +- 0.0%
  if (pt >= 200  && eta >= -2.750 && eta < -2.500) return 0.04009; // +- 0.0%
  if (pt >= 200  && eta >= -2.500 && eta < -2.250) return 0.02447; // +- 0.0%
  if (pt >= 200  && eta >= -2.250 && eta < -2.000) return 0.00839; // +- 0.0%
  if (pt >= 200  && eta >= 2.000 && eta < 2.250) return 0.01085; // +- 0.0%
  if (pt >= 200  && eta >= 2.250 && eta < 2.500) return 0.02865; // +- 0.0%
  if (pt >= 200  && eta >= 2.500 && eta < 2.750) return 0.03361; // +- 0.0%
  if (pt >= 200  && eta >= 2.750 && eta < 3.000) return 0.05753; // +- 0.0%
  return 0.;
}




std::vector<float> getPrefiringRates(float pt, float eta, int year, bool ispho) {
    // Dump maps using a ROOT file and the script from
    // /home/users/namin/2018/fourtop/all/FTAnalysis/analysis/checks/prefire/dump_scale_factors.py
    // year -- 2016 or 2017
    // ispho -- true if photon map otherwise jet map
    // returns central, up, down
    float prefiringRateSystUnc = 0.2;
    float rate = 0.;
    float error = 0.;
    if (ispho) {
        if (year == 2016) {
            rate = getPrefireInefficiency_singlephoton_2016(pt,eta);
            error = getPrefireInefficiencyError_singlephoton_2016(pt,eta);
        } else {
            rate = getPrefireInefficiency_singlephoton_2017(pt,eta);
            error = getPrefireInefficiencyError_singlephoton_2017(pt,eta);
        }
    } else {
        if (year == 2016) {
            rate = getPrefireInefficiency_singlejet_2016(pt,eta);
            error = getPrefireInefficiencyError_singlejet_2016(pt,eta);
        } else {
            rate = getPrefireInefficiency_singlejet_2017(pt,eta);
            error = getPrefireInefficiencyError_singlejet_2017(pt,eta);
        }
    }
    float rateup = min(max(rate+error, (float)(rate*(1.+prefiringRateSystUnc))),1.f);
    float ratedown = max(min(rate-error, (float)(rate*(1.-prefiringRateSystUnc))),0.f);
    // return std::make_tuple(rate, rateup, ratedown);
    return std::vector<float>({rate, rateup, ratedown});
}

std::tuple<float,float,float> getPrefireInfo(int year) {
    // CMSSW recipe details:
    //      https://twiki.cern.ch/twiki/bin/view/CMS/L1ECALPrefiringWeightRecipe
    // Implementation:
    //      https://github.com/lathomas/cmssw/blob/d6f39fee694f92a494ad7a4f4c7eb0433463335d/L1Prefiring/EventWeightProducer/plugins/L1ECALPrefiringWeightProducer.cc
    // Low efficiency for high eta EG objects in data before 2018 due to prefiring issue
    // This function looks at all jets/photons in an event and calculates a scale factor for simulation using
    //      https://github.com/lathomas/cmssw/blob/d6f39fee694f92a494ad7a4f4c7eb0433463335d/L1Prefiring/EventWeightProducer/files/L1PrefiringMaps_new.root
    // Protip: `std::tie(prefire_sf, prefire_sf_up, prefire_sf_down) = getPrefireInfo(gconf.year);`
    // More details in
    //     https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18
    // Sorry for the medium/big hardcoded lookup tables. It's slightly faster and smaller gzipped size than the root file.
    std::vector<float> nonprefiringprob = {1., 1., 1.};
    std::vector<LorentzVector> affectedphotons;
    for(unsigned int ipho = 0; ipho < cms3.photons_p4().size(); ipho++) {
        float pt = photons_p4()[ipho].pt();
        if (pt < 20.) continue;
        float eta = photons_p4()[ipho].eta();
        if (fabs(eta) < 2.) continue;
        if (fabs(eta) > 3.) continue;
        auto p4 = photons_p4()[ipho];
        affectedphotons.push_back(p4);
        // float cent, up, down;
        auto centupdown = getPrefiringRates(pt,eta,year,true);
        nonprefiringprob[0] *= 1.-centupdown[0];
        nonprefiringprob[1] *= 1.-centupdown[1];
        nonprefiringprob[2] *= 1.-centupdown[2];
    }
    for(unsigned int ijet = 0; ijet < cms3.pfjets_p4().size(); ijet++) {
        float jeteta = pfjets_p4()[ijet].eta();
        if (fabs(jeteta) < 2.) continue;
        if (fabs(jeteta) > 3.) continue;
        float jetpt = pfjets_p4()[ijet].pt(); //*pfjets_undoJEC()[ijet];
        if (jetpt < 20.) continue;

        float nonprefprob_pho_cent = 1.; // nonprefiring probability from overlapping photons
        float nonprefprob_pho_up = 1.;
        float nonprefprob_pho_down = 1.;
        auto jetp4 = pfjets_p4()[ijet];
        int noverlap = 0;
        for (const auto &phop4 : affectedphotons) {
            float phoeta = phop4.eta();
            if (DeltaR(jeteta,phoeta,phop4.phi(),jetp4.phi()) > 0.4) continue;
            float phopt = phop4.pt();
            auto centupdown = getPrefiringRates(phopt,phoeta,year,true); // prefiring probability of photon
            nonprefprob_pho_cent *= 1.-centupdown[0];
            nonprefprob_pho_up *= 1.-centupdown[1];
            nonprefprob_pho_down *= 1.-centupdown[2];
            noverlap += 1;
        }

        auto centupdown = getPrefiringRates(jetpt,jeteta,year,false); // prefiring probability of jet
        for (int i = 0; i < 3; i++) {
            if (noverlap == 0) {
                // If no overlapping photons, just multiply by jet nonprefiring rate
                nonprefiringprob[i] *= (1.-centupdown[i]);
            } else if (nonprefprob_pho_cent > 1.-centupdown[i]) {
                // If overlapping photons have a nonprefiring rate larger than the jet, replace those weights by the jet one
                nonprefiringprob[i] *= ((nonprefprob_pho_cent != 0.) ? (1.-centupdown[i])/nonprefprob_pho_cent : 0.);
            } else if (nonprefprob_pho_cent < 1.-centupdown[i]) {
                // If overlapping photons have a nonprefiring rate smaller than the jet, don't consider the jet in the event weight
                nonprefiringprob[i] *= 1.;
            }
        }
    }

    return std::make_tuple(nonprefiringprob[0], nonprefiringprob[1], nonprefiringprob[2]);
}
