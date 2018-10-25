#include "JetSelections.h"
#include "ElectronSelections.h"
#include "MuonSelections.h" 
#include "Math/VectorUtil.h"

using namespace tas;

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

float getPrefireInefficiency_singlejet(float pt, float eta) {
  if (pt >= 0 && pt < 10 && fabs(eta) >= 0.000 && fabs(eta) < 0.250) return 0.00002; // +- 57.7%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 0.250 && fabs(eta) < 0.500) return 0.00004; // +- 57.7%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 0.500 && fabs(eta) < 0.750) return 0.00001; // +- 100.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 0.750 && fabs(eta) < 1.000) return 0.00003; // +- 70.7%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 1.000 && fabs(eta) < 1.250) return 0.00003; // +- 70.7%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 1.250 && fabs(eta) < 1.500) return 0.00003; // +- 70.7%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00003; // +- 100.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00049; // +- 28.9%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00545; // +- 12.4%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00223; // +- 30.1%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00015; // +- 100.0%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00026; // +- 70.7%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00186; // +- 28.8%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00609; // +- 17.1%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00172; // +- 25.8%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00007; // +- 70.7%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00045; // +- 70.7%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00051; // +- 70.7%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00412; // +- 25.8%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01360; // +- 13.9%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00194; // +- 25.8%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00009; // +- 70.7%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 0.250 && fabs(eta) < 0.500) return 0.00013; // +- 100.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00028; // +- 100.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00118; // +- 50.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01499; // +- 13.9%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.04510; // +- 8.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01011; // +- 13.8%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00368; // +- 18.9%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00068; // +- 100.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00521; // +- 37.7%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.03828; // +- 14.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.12047; // +- 7.1%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.04659; // +- 10.4%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.05501; // +- 9.4%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.04762; // +- 97.6%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00278; // +- 70.6%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00281; // +- 70.6%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.05767; // +- 15.2%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.16435; // +- 8.4%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.08566; // +- 10.5%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.19619; // +- 6.8%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.20000; // +- 89.4%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.01955; // +- 37.4%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.07952; // +- 16.7%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.26136; // +- 8.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.21920; // +- 8.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.36822; // +- 5.8%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.25000; // +- 86.6%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.01282; // +- 49.7%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.08935; // +- 18.7%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.30960; // +- 8.3%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.39231; // +- 5.5%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.53778; // +- 4.4%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00806; // +- 99.6%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00602; // +- 99.7%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.08805; // +- 25.5%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.35912; // +- 9.9%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.47566; // +- 6.4%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.68093; // +- 4.3%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.02899; // +- 69.7%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.18919; // +- 24.1%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.50000; // +- 10.8%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.63846; // +- 6.6%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.78035; // +- 4.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.01724; // +- 99.1%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.01754; // +- 99.1%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.16949; // +- 28.8%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.48276; // +- 11.1%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.71818; // +- 6.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.80488; // +- 3.8%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.17241; // +- 40.7%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.32258; // +- 26.0%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.78049; // +- 8.3%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.90741; // +- 4.3%
  if (pt >= 300  && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.05882; // +- 97.0%
  if (pt >= 300  && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.14286; // +- 65.5%
  if (pt >= 300  && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.36364; // +- 39.9%
  if (pt >= 300  && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.90476; // +- 7.1%
  if (pt >= 300  && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.86047; // +- 6.1%
  return 0.;
}

float getPrefireInefficiencyError_singlejet(float pt, float eta) {
  if (pt >= 0 && pt < 10 && fabs(eta) >= 0.250 && fabs(eta) < 0.500) return 0.00002; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 0.500 && fabs(eta) < 0.750) return 0.00001; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 0.750 && fabs(eta) < 1.000) return 0.00002; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 1.000 && fabs(eta) < 1.250) return 0.00002; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 1.250 && fabs(eta) < 1.500) return 0.00002; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00003; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00014; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00067; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00067; // +- 0.0%
  if (pt >= 0 && pt < 10 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00015; // +- 0.0%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00018; // +- 0.0%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00054; // +- 0.0%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00104; // +- 0.0%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00044; // +- 0.0%
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00005; // +- 0.0%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00032; // +- 0.0%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00036; // +- 0.0%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00106; // +- 0.0%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00189; // +- 0.0%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00050; // +- 0.0%
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00006; // +- 0.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 0.250 && fabs(eta) < 0.500) return 0.00013; // +- 0.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00028; // +- 0.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00059; // +- 0.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00208; // +- 0.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00360; // +- 0.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00139; // +- 0.0%
  if (pt >= 30 && pt < 50 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00069; // +- 0.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00068; // +- 0.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00196; // +- 0.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00536; // +- 0.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00852; // +- 0.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00485; // +- 0.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00519; // +- 0.0%
  if (pt >= 50 && pt < 70 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.04647; // +- 0.0%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00196; // +- 0.0%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00199; // +- 0.0%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00874; // +- 0.0%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01383; // +- 0.0%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00899; // +- 0.0%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.01330; // +- 0.0%
  if (pt >= 70 && pt < 90 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.17889; // +- 0.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00732; // +- 0.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01328; // +- 0.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.02095; // +- 0.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01761; // +- 0.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.02123; // +- 0.0%
  if (pt >= 90 && pt < 110 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.21651; // +- 0.0%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00637; // +- 0.0%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01672; // +- 0.0%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.02572; // +- 0.0%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.02141; // +- 0.0%
  if (pt >= 110 && pt < 140 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.02350; // +- 0.0%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.00803; // +- 0.0%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.00601; // +- 0.0%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.02247; // +- 0.0%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.03566; // +- 0.0%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.03056; // +- 0.0%
  if (pt >= 140 && pt < 170 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.02908; // +- 0.0%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.02020; // +- 0.0%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.04553; // +- 0.0%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.05392; // +- 0.0%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.04214; // +- 0.0%
  if (pt >= 170 && pt < 200 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.03148; // +- 0.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 1.500 && fabs(eta) < 1.750) return 0.01709; // +- 0.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.01739; // +- 0.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.04884; // +- 0.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.05357; // +- 0.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.04289; // +- 0.0%
  if (pt >= 200 && pt < 250 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.03095; // +- 0.0%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.07014; // +- 0.0%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.08396; // +- 0.0%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.06464; // +- 0.0%
  if (pt >= 250 && pt < 300 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.03945; // +- 0.0%
  if (pt >= 300  && fabs(eta) >= 1.750 && fabs(eta) < 2.000) return 0.05707; // +- 0.0%
  if (pt >= 300  && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.09352; // +- 0.0%
  if (pt >= 300  && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.14504; // +- 0.0%
  if (pt >= 300  && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.06406; // +- 0.0%
  if (pt >= 300  && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.05284; // +- 0.0%
  return 0.;
}

std::tuple<float,float,int> getPrefireInfo() {
    // Low efficiency for high eta EG objects in data before 2018 due to prefiring issue
    // This function looks at all jets in an event and calculates a scale factor for simulation using
    //     https://ncsmith.web.cern.ch/ncsmith/PrefireEfficiencyMaps/Preliminary/Jet_L1IsoEG30eff_bxm1_looseJet_SingleMuon_Run2017F.pdf
    // The scale factor, error, and count of number of "affected" jets are returned as an std::tuple
    // Protip: `std::tie(prefire_sf, prefire_sferr, prefire_njets) = getPrefireInfo();`
    // More details in
    //     https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18
    //     https://github.com/nsmith-/PrefireAnalysis#jet-prefire-efficiencies
    // Sorry for the medium/big hardcoded lookup tables. The efficiency map is TEfficiency and I'm not going to read a ROOT file out of principle
    float sf = 1.;
    float sf_err = 0.; // error^2 for sf=(1-x1)(1-x2)...(1-xn) is (sf^2)*( (dx1/(1-x1))^2 + ... + (dx2/(1-x2))^2 )
    int naffected = 0;
    for(unsigned int ijet = 0; ijet < cms3.pfjets_p4().size(); ijet++) {
        float jeteta = pfjets_p4()[ijet].eta();
        // skip for performance; these jets have 0 or negligible inefficiency
        if (fabs(jeteta) < 1.5) continue;
        // calculate raw EM jet pT
        float jetempt = pfjets_p4()[ijet].pt()*(pfjets_chargedEmE()[ijet]+pfjets_neutralEmE()[ijet])/pfjets_p4()[ijet].E();
        float ineff = getPrefireInefficiency_singlejet(jetempt, jeteta);
        float ineff_err = getPrefireInefficiencyError_singlejet(jetempt, jeteta);
        sf *= 1.-ineff;
        sf_err += pow(ineff_err/(1-ineff),2);
        naffected += (ineff > 1e-6);
    }
    sf_err = sf*pow(sf_err,0.5);
    return {sf, sf_err, naffected};
}
