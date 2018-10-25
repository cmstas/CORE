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

float getPrefireInefficiency_singlejet_2016(float pt, float eta) {
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01168; // +- 6.7%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01192; // +- 6.6%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01274; // +- 6.5%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.01397; // +- 5.6%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01109; // +- 9.3%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01303; // +- 7.2%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01324; // +- 7.4%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01174; // +- 8.1%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.01162; // +- 8.0%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01320; // +- 10.6%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01172; // +- 8.2%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01542; // +- 7.4%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01427; // +- 8.4%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.01071; // +- 9.8%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01266; // +- 12.6%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01255; // +- 8.8%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01821; // +- 7.7%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01570; // +- 9.1%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.01190; // +- 11.3%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01213; // +- 14.7%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01376; // +- 10.0%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.02407; // +- 8.1%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.02423; // +- 8.8%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.02252; // +- 10.3%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00887; // +- 21.2%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01463; // +- 10.3%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.03302; // +- 7.5%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.03006; // +- 8.6%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.02823; // +- 10.3%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01101; // +- 22.2%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01786; // +- 10.8%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.04995; // +- 7.1%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.05086; // +- 7.9%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.05818; // +- 8.6%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01060; // +- 27.6%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01830; // +- 11.7%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.06205; // +- 7.2%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.08718; // +- 6.7%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.11075; // +- 7.3%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01129; // +- 33.1%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.02524; // +- 10.7%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.07743; // +- 7.0%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.14886; // +- 5.7%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.17355; // +- 6.6%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.02104; // +- 27.4%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.02197; // +- 12.3%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.09868; // +- 6.6%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.20642; // +- 5.3%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.28276; // +- 5.9%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01034; // +- 49.7%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.03256; // +- 12.0%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.12474; // +- 6.9%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.30693; // +- 5.0%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.34579; // +- 6.6%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01807; // +- 57.2%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.02860; // +- 17.2%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.15698; // +- 7.9%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.35200; // +- 6.1%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.50000; // +- 7.3%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.02857; // +- 69.7%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.04407; // +- 15.7%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.15918; // +- 10.4%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.41104; // +- 6.6%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.50495; // +- 9.9%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.04469; // +- 34.6%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.20290; // +- 23.9%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.64706; // +- 12.7%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.91667; // +- 8.7%
  if (pt >= 574 && pt < 830 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.33333; // +- 47.1%
  if (pt >= 574 && pt < 830 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 1.00000; // +- 0.0%
  return 0.;
}

float getPrefireInefficiencyError_singlejet_2016(float pt, float eta) {
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00079; // +- 0.0%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00079; // +- 0.0%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00083; // +- 0.0%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00079; // +- 0.0%
  if (pt >= 30 && pt < 36 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00103; // +- 0.0%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00094; // +- 0.0%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00098; // +- 0.0%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00095; // +- 0.0%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00092; // +- 0.0%
  if (pt >= 36 && pt < 43 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00140; // +- 0.0%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00096; // +- 0.0%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00114; // +- 0.0%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00119; // +- 0.0%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00106; // +- 0.0%
  if (pt >= 43 && pt < 52 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00160; // +- 0.0%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00111; // +- 0.0%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00141; // +- 0.0%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00142; // +- 0.0%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00135; // +- 0.0%
  if (pt >= 52 && pt < 63 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00178; // +- 0.0%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00138; // +- 0.0%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00195; // +- 0.0%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00214; // +- 0.0%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00231; // +- 0.0%
  if (pt >= 63 && pt < 75 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00188; // +- 0.0%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00151; // +- 0.0%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00247; // +- 0.0%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00259; // +- 0.0%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00290; // +- 0.0%
  if (pt >= 75 && pt < 91 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00245; // +- 0.0%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00193; // +- 0.0%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00357; // +- 0.0%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00399; // +- 0.0%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00501; // +- 0.0%
  if (pt >= 91 && pt < 109 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00293; // +- 0.0%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00214; // +- 0.0%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00444; // +- 0.0%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00583; // +- 0.0%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.00803; // +- 0.0%
  if (pt >= 109 && pt < 131 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00374; // +- 0.0%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00270; // +- 0.0%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00544; // +- 0.0%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.00848; // +- 0.0%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.01148; // +- 0.0%
  if (pt >= 131 && pt < 158 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00577; // +- 0.0%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00270; // +- 0.0%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00647; // +- 0.0%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01093; // +- 0.0%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.01673; // +- 0.0%
  if (pt >= 158 && pt < 190 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.00514; // +- 0.0%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00391; // +- 0.0%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.00865; // +- 0.0%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.01530; // +- 0.0%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.02299; // +- 0.0%
  if (pt >= 190 && pt < 228 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01034; // +- 0.0%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00491; // +- 0.0%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01240; // +- 0.0%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.02136; // +- 0.0%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.03666; // +- 0.0%
  if (pt >= 228 && pt < 274 && fabs(eta) >= 3.000 && fabs(eta) < 3.500) return 0.01991; // +- 0.0%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.00690; // +- 0.0%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.01653; // +- 0.0%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.02725; // +- 0.0%
  if (pt >= 274 && pt < 397 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.04975; // +- 0.0%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.000 && fabs(eta) < 2.250) return 0.01544; // +- 0.0%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.04841; // +- 0.0%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.500 && fabs(eta) < 2.750) return 0.08196; // +- 0.0%
  if (pt >= 397 && pt < 574 && fabs(eta) >= 2.750 && fabs(eta) < 3.000) return 0.07979; // +- 0.0%
  if (pt >= 574 && pt < 830 && fabs(eta) >= 2.250 && fabs(eta) < 2.500) return 0.15713; // +- 0.0%
  return 0.;
}

float getPrefireInefficiency_singlejet_2017(float pt, float eta) {
  if (pt >= 40 && pt < 50 && eta >= -3.500 && eta < -3.250) return 0.00616; // +- 11.2%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.00865; // +- 10.5%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.01134; // +- 12.6%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.00940; // +- 13.9%
  if (pt >= 40 && pt < 50 && eta >= 1.750 && eta < 2.000) return 0.00531; // +- 5.9%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.00710; // +- 12.5%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.01452; // +- 13.2%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.01031; // +- 12.1%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.00865; // +- 9.3%
  if (pt >= 40 && pt < 50 && eta >= 3.000 && eta < 3.250) return 0.00524; // +- 6.2%
  if (pt >= 50 && pt < 60 && eta >= -3.500 && eta < -3.250) return 0.00511; // +- 5.9%
  if (pt >= 50 && pt < 60 && eta >= -3.000 && eta < -2.750) return 0.01402; // +- 15.7%
  if (pt >= 50 && pt < 60 && eta >= -2.750 && eta < -2.500) return 0.01121; // +- 18.0%
  if (pt >= 50 && pt < 60 && eta >= -2.500 && eta < -2.250) return 0.02075; // +- 15.5%
  if (pt >= 50 && pt < 60 && eta >= 2.000 && eta < 2.250) return 0.00615; // +- 13.7%
  if (pt >= 50 && pt < 60 && eta >= 2.250 && eta < 2.500) return 0.02192; // +- 15.4%
  if (pt >= 50 && pt < 60 && eta >= 2.500 && eta < 2.750) return 0.01959; // +- 14.7%
  if (pt >= 50 && pt < 60 && eta >= 2.750 && eta < 3.000) return 0.01626; // +- 14.4%
  if (pt >= 60 && pt < 70 && eta >= -3.250 && eta < -3.000) return 0.00946; // +- 28.0%
  if (pt >= 60 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.02612; // +- 18.2%
  if (pt >= 60 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.02831; // +- 19.1%
  if (pt >= 60 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.02661; // +- 17.5%
  if (pt >= 60 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.00636; // +- 18.9%
  if (pt >= 60 && pt < 70 && eta >= -2.000 && eta < -1.750) return 0.00540; // +- 11.1%
  if (pt >= 60 && pt < 70 && eta >= 1.750 && eta < 2.000) return 0.00746; // +- 19.1%
  if (pt >= 60 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.01116; // +- 22.3%
  if (pt >= 60 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.05349; // +- 13.3%
  if (pt >= 60 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.02523; // +- 17.7%
  if (pt >= 60 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.03247; // +- 15.6%
  if (pt >= 60 && pt < 70 && eta >= 3.000 && eta < 3.250) return 0.00789; // +- 27.0%
  if (pt >= 70 && pt < 80 && eta >= -3.250 && eta < -3.000) return 0.00787; // +- 34.8%
  if (pt >= 70 && pt < 80 && eta >= -3.000 && eta < -2.750) return 0.05905; // +- 16.7%
  if (pt >= 70 && pt < 80 && eta >= -2.750 && eta < -2.500) return 0.04255; // +- 18.8%
  if (pt >= 70 && pt < 80 && eta >= -2.500 && eta < -2.250) return 0.03377; // +- 21.4%
  if (pt >= 70 && pt < 80 && eta >= -2.250 && eta < -2.000) return 0.00621; // +- 22.1%
  if (pt >= 70 && pt < 80 && eta >= 1.750 && eta < 2.000) return 0.00922; // +- 25.5%
  if (pt >= 70 && pt < 80 && eta >= 2.000 && eta < 2.250) return 0.01294; // +- 27.6%
  if (pt >= 70 && pt < 80 && eta >= 2.250 && eta < 2.500) return 0.06111; // +- 16.2%
  if (pt >= 70 && pt < 80 && eta >= 2.500 && eta < 2.750) return 0.04348; // +- 18.4%
  if (pt >= 70 && pt < 80 && eta >= 2.750 && eta < 3.000) return 0.05009; // +- 17.5%
  if (pt >= 70 && pt < 80 && eta >= 3.000 && eta < 3.250) return 0.00510; // +- 10.0%
  if (pt >= 70 && pt < 80 && eta >= 3.250 && eta < 3.500) return 0.00847; // +- 36.9%
  if (pt >= 80 && pt < 90 && eta >= -3.000 && eta < -2.750) return 0.08157; // +- 17.9%
  if (pt >= 80 && pt < 90 && eta >= -2.750 && eta < -2.500) return 0.06540; // +- 19.0%
  if (pt >= 80 && pt < 90 && eta >= -2.500 && eta < -2.250) return 0.06497; // +- 19.4%
  if (pt >= 80 && pt < 90 && eta >= -2.250 && eta < -2.000) return 0.01515; // +- 33.2%
  if (pt >= 80 && pt < 90 && eta >= -2.000 && eta < -1.750) return 0.00771; // +- 29.6%
  if (pt >= 80 && pt < 90 && eta >= 1.750 && eta < 2.000) return 0.00821; // +- 31.2%
  if (pt >= 80 && pt < 90 && eta >= 2.000 && eta < 2.250) return 0.01896; // +- 30.1%
  if (pt >= 80 && pt < 90 && eta >= 2.250 && eta < 2.500) return 0.09065; // +- 16.4%
  if (pt >= 80 && pt < 90 && eta >= 2.500 && eta < 2.750) return 0.08683; // +- 16.7%
  if (pt >= 80 && pt < 90 && eta >= 2.750 && eta < 3.000) return 0.09249; // +- 16.4%
  if (pt >= 80 && pt < 90 && eta >= 3.000 && eta < 3.250) return 0.00810; // +- 43.7%
  if (pt >= 90 && pt < 100 && eta >= -3.500 && eta < -3.250) return 0.01333; // +- 55.7%
  if (pt >= 90 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.12565; // +- 18.8%
  if (pt >= 90 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.11345; // +- 17.8%
  if (pt >= 90 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.05213; // +- 28.0%
  if (pt >= 90 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.01805; // +- 37.8%
  if (pt >= 90 && pt < 100 && eta >= -2.000 && eta < -1.750) return 0.00581; // +- 26.4%
  if (pt >= 90 && pt < 100 && eta >= 1.750 && eta < 2.000) return 0.00915; // +- 38.8%
  if (pt >= 90 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.00800; // +- 43.2%
  if (pt >= 90 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.11161; // +- 18.5%
  if (pt >= 90 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.09562; // +- 18.9%
  if (pt >= 90 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.14091; // +- 16.4%
  if (pt >= 90 && pt < 100 && eta >= 3.000 && eta < 3.250) return 0.00654; // +- 48.4%
  if (pt >= 100 && pt < 125 && eta >= -3.250 && eta < -3.000) return 0.00546; // +- 29.1%
  if (pt >= 100 && pt < 125 && eta >= -3.000 && eta < -2.750) return 0.30769; // +- 9.8%
  if (pt >= 100 && pt < 125 && eta >= -2.750 && eta < -2.500) return 0.16608; // +- 13.2%
  if (pt >= 100 && pt < 125 && eta >= -2.500 && eta < -2.250) return 0.10667; // +- 16.4%
  if (pt >= 100 && pt < 125 && eta >= -2.250 && eta < -2.000) return 0.02116; // +- 30.6%
  if (pt >= 100 && pt < 125 && eta >= 1.750 && eta < 2.000) return 0.00682; // +- 29.8%
  if (pt >= 100 && pt < 125 && eta >= 2.000 && eta < 2.250) return 0.04776; // +- 23.1%
  if (pt >= 100 && pt < 125 && eta >= 2.250 && eta < 2.500) return 0.12579; // +- 14.5%
  if (pt >= 100 && pt < 125 && eta >= 2.500 && eta < 2.750) return 0.19016; // +- 11.7%
  if (pt >= 100 && pt < 125 && eta >= 2.750 && eta < 3.000) return 0.33203; // +- 8.8%
  if (pt >= 100 && pt < 125 && eta >= 3.000 && eta < 3.250) return 0.00543; // +- 28.3%
  if (pt >= 100 && pt < 125 && eta >= 3.250 && eta < 3.500) return 0.00649; // +- 47.9%
  if (pt >= 125 && pt < 150 && eta >= -3.500 && eta < -3.250) return 0.03390; // +- 64.3%
  if (pt >= 125 && pt < 150 && eta >= -3.250 && eta < -3.000) return 0.01587; // +- 82.3%
  if (pt >= 125 && pt < 150 && eta >= -3.000 && eta < -2.750) return 0.46429; // +- 11.7%
  if (pt >= 125 && pt < 150 && eta >= -2.750 && eta < -2.500) return 0.28462; // +- 13.8%
  if (pt >= 125 && pt < 150 && eta >= -2.500 && eta < -2.250) return 0.14388; // +- 20.4%
  if (pt >= 125 && pt < 150 && eta >= -2.250 && eta < -2.000) return 0.03593; // +- 37.3%
  if (pt >= 125 && pt < 150 && eta >= -2.000 && eta < -1.750) return 0.01579; // +- 47.5%
  if (pt >= 125 && pt < 150 && eta >= 1.750 && eta < 2.000) return 0.00971; // +- 49.1%
  if (pt >= 125 && pt < 150 && eta >= 2.000 && eta < 2.250) return 0.06211; // +- 29.4%
  if (pt >= 125 && pt < 150 && eta >= 2.250 && eta < 2.500) return 0.23129; // +- 14.9%
  if (pt >= 125 && pt < 150 && eta >= 2.500 && eta < 2.750) return 0.26087; // +- 14.2%
  if (pt >= 125 && pt < 150 && eta >= 2.750 && eta < 3.000) return 0.50420; // +- 9.1%
  if (pt >= 125 && pt < 150 && eta >= 3.000 && eta < 3.250) return 0.02740; // +- 63.2%
  if (pt >= 125 && pt < 150 && eta >= 3.250 && eta < 3.500) return 0.02174; // +- 87.0%
  if (pt >= 150 && pt < 200 && eta >= -3.250 && eta < -3.000) return 0.02439; // +- 88.3%
  if (pt >= 150 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.67532; // +- 7.9%
  if (pt >= 150 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.31111; // +- 15.6%
  if (pt >= 150 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.16071; // +- 21.3%
  if (pt >= 150 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.01527; // +- 57.7%
  if (pt >= 150 && pt < 200 && eta >= 1.750 && eta < 2.000) return 0.01987; // +- 49.6%
  if (pt >= 150 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.03788; // +- 41.0%
  if (pt >= 150 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.22727; // +- 19.5%
  if (pt >= 150 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.52066; // +- 8.7%
  if (pt >= 150 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.62857; // +- 9.2%
  if (pt >= 150 && pt < 200 && eta >= 3.000 && eta < 3.250) return 0.01754; // +- 84.0%
  if (pt >= 150 && pt < 200 && eta >= 3.250 && eta < 3.500) return 0.02632; // +- 89.0%
  if (pt >= 200 && pt < 300 && eta >= -3.000 && eta < -2.750) return 0.81818; // +- 10.2%
  if (pt >= 200 && pt < 300 && eta >= -2.750 && eta < -2.500) return 0.40741; // +- 23.2%
  if (pt >= 200 && pt < 300 && eta >= -2.500 && eta < -2.250) return 0.17073; // +- 34.0%
  if (pt >= 200 && pt < 300 && eta >= -2.250 && eta < -2.000) return 0.03509; // +- 64.5%
  if (pt >= 200 && pt < 300 && eta >= 1.750 && eta < 2.000) return 0.04348; // +- 53.3%
  if (pt >= 200 && pt < 300 && eta >= 2.000 && eta < 2.250) return 0.11290; // +- 34.9%
  if (pt >= 200 && pt < 300 && eta >= 2.250 && eta < 2.500) return 0.39394; // +- 21.5%
  if (pt >= 200 && pt < 300 && eta >= 2.500 && eta < 2.750) return 0.63265; // +- 10.9%
  if (pt >= 200 && pt < 300 && eta >= 2.750 && eta < 3.000) return 0.78788; // +- 9.1%
  if (pt >= 200 && pt < 300 && eta >= 3.000 && eta < 3.250) return 0.06667; // +- 93.2%
  if (pt >= 300  && eta >= -3.000 && eta < -2.750) return 1.00000; // +- 2.7%
  if (pt >= 300  && eta >= -2.750 && eta < -2.500) return 0.85714; // +- 15.7%
  if (pt >= 300  && eta >= -2.500 && eta < -2.250) return 0.28571; // +- 59.4%
  if (pt >= 300  && eta >= -2.250 && eta < -2.000) return 0.11111; // +- 92.4%
  if (pt >= 300  && eta >= 2.000 && eta < 2.250) return 0.11111; // +- 92.4%
  if (pt >= 300  && eta >= 2.250 && eta < 2.500) return 0.58333; // +- 24.4%
  if (pt >= 300  && eta >= 2.500 && eta < 2.750) return 0.60000; // +- 36.6%
  if (pt >= 300  && eta >= 2.750 && eta < 3.000) return 0.85714; // +- 15.7%
  return 0.;
}

float getPrefireInefficiencyError_singlejet_2017(float pt, float eta) {
  if (pt >= 40 && pt < 50 && eta >= -3.500 && eta < -3.250) return 0.00069; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.250 && eta < -3.000) return 0.00137; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -3.000 && eta < -2.750) return 0.00091; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.750 && eta < -2.500) return 0.00143; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.250) return 0.00131; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.250 && eta < -2.000) return 0.00131; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= -2.000 && eta < -1.750) return 0.00122; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 1.750 && eta < 2.000) return 0.00031; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.250) return 0.00089; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.250 && eta < 2.500) return 0.00192; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.500 && eta < 2.750) return 0.00125; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 2.750 && eta < 3.000) return 0.00080; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 3.000 && eta < 3.250) return 0.00033; // +- 0.0%
  if (pt >= 40 && pt < 50 && eta >= 3.250 && eta < 3.500) return 0.00150; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= -3.500 && eta < -3.250) return 0.00030; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= -3.250 && eta < -3.000) return 0.00205; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= -3.000 && eta < -2.750) return 0.00220; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= -2.750 && eta < -2.500) return 0.00202; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= -2.500 && eta < -2.250) return 0.00322; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= -2.250 && eta < -2.000) return 0.00159; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= -2.000 && eta < -1.750) return 0.00151; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= 1.750 && eta < 2.000) return 0.00157; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= 2.000 && eta < 2.250) return 0.00084; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= 2.250 && eta < 2.500) return 0.00338; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= 2.500 && eta < 2.750) return 0.00288; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= 2.750 && eta < 3.000) return 0.00234; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= 3.000 && eta < 3.250) return 0.00204; // +- 0.0%
  if (pt >= 50 && pt < 60 && eta >= 3.250 && eta < 3.500) return 0.00183; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= -3.500 && eta < -3.250) return 0.00205; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= -3.250 && eta < -3.000) return 0.00265; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= -3.000 && eta < -2.750) return 0.00474; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= -2.750 && eta < -2.500) return 0.00541; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= -2.500 && eta < -2.250) return 0.00465; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= -2.250 && eta < -2.000) return 0.00120; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= -2.000 && eta < -1.750) return 0.00060; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= 1.750 && eta < 2.000) return 0.00143; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= 2.000 && eta < 2.250) return 0.00249; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= 2.250 && eta < 2.500) return 0.00710; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= 2.500 && eta < 2.750) return 0.00447; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= 2.750 && eta < 3.000) return 0.00505; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= 3.000 && eta < 3.250) return 0.00213; // +- 0.0%
  if (pt >= 60 && pt < 70 && eta >= 3.250 && eta < 3.500) return 0.00292; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= -3.500 && eta < -3.250) return 0.00287; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= -3.250 && eta < -3.000) return 0.00274; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= -3.000 && eta < -2.750) return 0.00987; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= -2.750 && eta < -2.500) return 0.00801; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= -2.500 && eta < -2.250) return 0.00724; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= -2.250 && eta < -2.000) return 0.00137; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= -2.000 && eta < -1.750) return 0.00254; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= 1.750 && eta < 2.000) return 0.00235; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= 2.000 && eta < 2.250) return 0.00357; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= 2.250 && eta < 2.500) return 0.00990; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= 2.500 && eta < 2.750) return 0.00802; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= 2.750 && eta < 3.000) return 0.00878; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= 3.000 && eta < 3.250) return 0.00051; // +- 0.0%
  if (pt >= 70 && pt < 80 && eta >= 3.250 && eta < 3.500) return 0.00313; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= -3.500 && eta < -3.250) return 0.00414; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= -3.250 && eta < -3.000) return 0.00620; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= -3.000 && eta < -2.750) return 0.01462; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= -2.750 && eta < -2.500) return 0.01243; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= -2.500 && eta < -2.250) return 0.01262; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= -2.250 && eta < -2.000) return 0.00504; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= -2.000 && eta < -1.750) return 0.00228; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= 1.750 && eta < 2.000) return 0.00256; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= 2.000 && eta < 2.250) return 0.00571; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= 2.250 && eta < 2.500) return 0.01489; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= 2.500 && eta < 2.750) return 0.01451; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= 2.750 && eta < 3.000) return 0.01519; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= 3.000 && eta < 3.250) return 0.00354; // +- 0.0%
  if (pt >= 80 && pt < 90 && eta >= 3.250 && eta < 3.500) return 0.00468; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= -3.500 && eta < -3.250) return 0.00742; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= -3.250 && eta < -3.000) return 0.00760; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= -3.000 && eta < -2.750) return 0.02357; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= -2.750 && eta < -2.500) return 0.02016; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= -2.500 && eta < -2.250) return 0.01459; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= -2.250 && eta < -2.000) return 0.00682; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= -2.000 && eta < -1.750) return 0.00154; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= 1.750 && eta < 2.000) return 0.00355; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= 2.000 && eta < 2.250) return 0.00346; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= 2.250 && eta < 2.500) return 0.02062; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= 2.500 && eta < 2.750) return 0.01812; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= 2.750 && eta < 3.000) return 0.02310; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= 3.000 && eta < 3.250) return 0.00317; // +- 0.0%
  if (pt >= 90 && pt < 100 && eta >= 3.250 && eta < 3.500) return 0.00816; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= -3.500 && eta < -3.250) return 0.00522; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= -3.250 && eta < -3.000) return 0.00159; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= -3.000 && eta < -2.750) return 0.03003; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= -2.750 && eta < -2.500) return 0.02185; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= -2.500 && eta < -2.250) return 0.01745; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= -2.250 && eta < -2.000) return 0.00649; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= -2.000 && eta < -1.750) return 0.00334; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= 1.750 && eta < 2.000) return 0.00203; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= 2.000 && eta < 2.250) return 0.01105; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= 2.250 && eta < 2.500) return 0.01827; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= 2.500 && eta < 2.750) return 0.02224; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= 2.750 && eta < 3.000) return 0.02932; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= 3.000 && eta < 3.250) return 0.00154; // +- 0.0%
  if (pt >= 100 && pt < 125 && eta >= 3.250 && eta < 3.500) return 0.00311; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= -3.500 && eta < -3.250) return 0.02181; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= -3.250 && eta < -3.000) return 0.01307; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= -3.000 && eta < -2.750) return 0.05437; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= -2.750 && eta < -2.500) return 0.03936; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= -2.500 && eta < -2.250) return 0.02933; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= -2.250 && eta < -2.000) return 0.01340; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= -2.000 && eta < -1.750) return 0.00749; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= 1.750 && eta < 2.000) return 0.00477; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= 2.000 && eta < 2.250) return 0.01829; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= 2.250 && eta < 2.500) return 0.03451; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= 2.500 && eta < 2.750) return 0.03714; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= 2.750 && eta < 3.000) return 0.04583; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= 3.000 && eta < 3.250) return 0.01732; // +- 0.0%
  if (pt >= 125 && pt < 150 && eta >= 3.250 && eta < 3.500) return 0.01892; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -3.500 && eta < -3.250) return 0.02898; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -3.250 && eta < -3.000) return 0.02154; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -3.000 && eta < -2.750) return 0.05357; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.750 && eta < -2.500) return 0.04858; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.500 && eta < -2.250) return 0.03426; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.250 && eta < -2.000) return 0.00881; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= -2.000 && eta < -1.750) return 0.00604; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 1.750 && eta < 2.000) return 0.00985; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.000 && eta < 2.250) return 0.01552; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.250 && eta < 2.500) return 0.04432; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.500 && eta < 2.750) return 0.04543; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 2.750 && eta < 3.000) return 0.05791; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 3.000 && eta < 3.250) return 0.01474; // +- 0.0%
  if (pt >= 150 && pt < 200 && eta >= 3.250 && eta < 3.500) return 0.02343; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -3.500 && eta < -3.250) return 0.10476; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -3.250 && eta < -3.000) return 0.04647; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -3.000 && eta < -2.750) return 0.08310; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.750 && eta < -2.500) return 0.09437; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.500 && eta < -2.250) return 0.05807; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.250 && eta < -2.000) return 0.02263; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= -2.000 && eta < -1.750) return 0.01379; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 1.750 && eta < 2.000) return 0.02316; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.000 && eta < 2.250) return 0.03940; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.250 && eta < 2.500) return 0.08486; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.500 && eta < 2.750) return 0.06906; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 2.750 && eta < 3.000) return 0.07177; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 3.000 && eta < 3.250) return 0.06211; // +- 0.0%
  if (pt >= 200 && pt < 300 && eta >= 3.250 && eta < 3.500) return 0.06441; // +- 0.0%
  if (pt >= 300  && eta >= -3.250 && eta < -3.000) return 0.21651; // +- 0.0%
  if (pt >= 300  && eta >= -3.000 && eta < -2.750) return 0.02666; // +- 0.0%
  if (pt >= 300  && eta >= -2.750 && eta < -2.500) return 0.13416; // +- 0.0%
  if (pt >= 300  && eta >= -2.500 && eta < -2.250) return 0.16984; // +- 0.0%
  if (pt >= 300  && eta >= -2.250 && eta < -2.000) return 0.10266; // +- 0.0%
  if (pt >= 300  && eta >= -2.000 && eta < -1.750) return 0.07391; // +- 0.0%
  if (pt >= 300  && eta >= 1.750 && eta < 2.000) return 0.04441; // +- 0.0%
  if (pt >= 300  && eta >= 2.000 && eta < 2.250) return 0.10266; // +- 0.0%
  if (pt >= 300  && eta >= 2.250 && eta < 2.500) return 0.14256; // +- 0.0%
  if (pt >= 300  && eta >= 2.500 && eta < 2.750) return 0.21953; // +- 0.0%
  if (pt >= 300  && eta >= 2.750 && eta < 3.000) return 0.13416; // +- 0.0%
  if (pt >= 300  && eta >= 3.000 && eta < 3.250) return 0.21651; // +- 0.0%
  return 0.;
}

std::tuple<float,float,int> getPrefireInfo(int year) {
    // Low efficiency for high eta EG objects in data before 2018 due to prefiring issue
    // This function looks at all jets in an event and calculates a scale factor for simulation using
    //     if year == 2016: https://ncsmith.web.cern.ch/ncsmith/PrefireEfficiencyMaps/Preliminary/Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.pdf (and .root)
    //     if year == 2017: Jet map 2017BtoF from https://lathomas.web.cern.ch/lathomas/TSGStuff/L1Prefiring/PrefiringMaps/
    // The scale factor, error, and count of number of "affected" jets are returned as an std::tuple
    // Protip: `std::tie(prefire_sf, prefire_sferr, prefire_njets) = getPrefireInfo(gconf.year);`
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
        if (fabs(jeteta) < 1.75) continue;
        // calculate raw EM jet pT
        // float jetempt = pfjets_p4()[ijet].pt()*(pfjets_chargedEmE()[ijet]+pfjets_neutralEmE()[ijet])/pfjets_p4()[ijet].E();
        // new maps use jet pT
        float jetpt = pfjets_p4()[ijet].pt();
        float ineff = year == 2016 ? getPrefireInefficiency_singlejet_2016(jetpt,jeteta) : getPrefireInefficiency_singlejet_2017(jetpt,jeteta);
        float ineff_err = year == 2016 ? getPrefireInefficiencyError_singlejet_2016(jetpt,jeteta) : getPrefireInefficiencyError_singlejet_2017(jetpt,jeteta);
        sf *= 1.-ineff;
        sf_err += pow(ineff_err/(1-ineff),2);
        naffected += (ineff > 1e-6);
    }
    sf_err = sf*pow(sf_err,0.5);
    return {sf, sf_err, naffected};
}
