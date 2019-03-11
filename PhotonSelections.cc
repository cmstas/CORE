#include <stdexcept>

#include "PhotonSelections.h"
#include "IsolationTools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

using namespace tas;

bool isLoosePhoton(unsigned int phIdx, analysis_t analysis, int version){
  // This is based on POG LOOSE, with an additional sieie sideband
  if (analysis == HAD){
    if (version == 1){
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.015) return false; // looser than POG (0.01) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.0559) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.035) return false; // larger than POG (0.0321) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.049) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else return false;
    }
    if (version == 2){
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.015) return false; // looser than POG (0.01) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.048) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.035) return false; // larger than POG (0.0321) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.069) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else return false;
    } 
    if (version == 3){ // Based on 50ns Spring15 IDs
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.015) return false; // looser than POG (0.01) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.05) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.035) return false; // larger than POG (0.0321) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.05) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else return false;
    } 
  }
  else {
    std::cerr << "CORE::PhotonSelection -- Photon selection not implemented for analysis " << analysis << std::endl;
    std::cerr << "CORE::PhotonSelection -- I will return true." << std::endl;
  }

  return true;
}

bool isTightPhoton(unsigned int phIdx, analysis_t analysis, int version){
  // This should be the same as POG LOOSE: 
  //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
  if (analysis == HAD){
    if (version == 1){
      if (!isLoosePhoton(phIdx, HAD)) return false;
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.010) return false; 
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0321) return false; 
        return true;
      } 
      else return false;
    }
    if (version == 2){
      if (!isLoosePhoton(phIdx, HAD, 2)) return false;
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0106) return false; 
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0266) return false; 
        return true;
      } 
      else return false;
    } 
    if (version == 3){
      if (!isLoosePhoton(phIdx, HAD, 3)) return false;
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0103) return false; 
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0277) return false; 
        return true;
      } 
      else return false;
    } 
  }
  else {
    std::cerr << "CORE::PhotonSelection -- Photon selection not implemented for analysis " << analysis << std::endl;
    std::cerr << "CORE::PhotonSelection -- I will return true." << std::endl;
  } 

  return true;
}

bool photonID(unsigned int phIdx, id_level_t id_level){
    
  // analysis_t analysis = whichAnalysis(id_level);

  switch (id_level){

	///////////////////
	// MET Templates //
	///////////////////

  case(ZMET_photon_v5):
    if( !isTemplatePhoton( phIdx ) ) return false;
	if( !isLoosePhotonPOG_Fall17V2( phIdx ) ) return false;
        // add trigger emulation cuts here, compared to v3
	if( !passTriggerEmu( phIdx ) ) return false;
	// Following cuts done in analysis code.
	// match to pfjet w/ neutral EM fraction > 70%
	// reject photons within electron with pT > 10 within cone of dR < 0.2
	// Reject photons aligned to MET within 0.14 radians in phi
	else return true;

      break;
  case(ZMET_photon_v4):

	if( !isTemplatePhoton( phIdx ) ) return false;
	if( !isLoosePhoton_Spring15_25ns( phIdx ) ) return false;
        // add trigger emulation cuts here, compared to v3
	if( !passTriggerEmu( phIdx ) ) return false;
	// Following cuts done in analysis code.
	// match to pfjet w/ neutral EM fraction > 70%
	// reject photons within electron with pT > 10 within cone of dR < 0.2
	// Reject photons aligned to MET within 0.14 radians in phi
	else return true;
	break;
  
  case(ZMET_photon_v3):

	if( !isTemplatePhoton( phIdx ) ) return false;
	if( !isLoosePhoton_Spring15_25ns( phIdx ) ) return false;
	// Following cuts done in analysis code.
	// match to pfjet w/ neutral EM fraction > 70%
	// reject photons within electron with pT > 10 within cone of dR < 0.2
	// Reject photons aligned to MET within 0.14 radians in phi
	else return true;
	break;
  
  case(ZMET_photon_v2):

	if( !isTemplatePhoton( phIdx ) ) return false;
	if( !isLoosePhoton_Spring15_50ns( phIdx ) ) return false;
	// Following cuts done in analysis code.
	// match to pfjet w/ neutral EM fraction > 70%
	// reject photons within electron with pT > 10 within cone of dR < 0.2
	// Reject photons aligned to MET within 0.14 radians in phi
	else return true;
	break;
  
  case(ZMET_photon_v1):

	if( !isTemplatePhoton( phIdx ) ) return false;
	// Following cuts done in analysis code.
	// match to pfjet w/ neutral EM fraction > 70%
	// reject photons within electron with pT > 10 within cone of dR < 0.2
	// Reject photons aligned to MET within 0.14 radians in phi
	else return true;
	break;

  case(STOP_loose_v4):
    if( gconf.year == 2016 ){
      if( !isLoosePhotonPOG_Spring16( phIdx ) ) return false;
    }
    else if( gconf.year >= 2017 ){
      if( !isLoosePhotonPOG_Fall17V2( phIdx ) ) return false;
    }
    return true;
    break;

  case(STOP_medium_v4):
    if( gconf.year == 2016 ){
      if( !isMediumPhotonPOG_Spring16( phIdx ) ) return false;
    }
    else if( gconf.year >= 2017 ){
      if( !isMediumPhotonPOG_Fall17V2( phIdx ) ) return false;
    }
    return true;
    break;

  case(STOP_tight_v4):
    if( gconf.year == 2016 ){
      if( !isTightPhotonPOG_Spring16( phIdx ) ) return false;
    }
    else if( gconf.year >= 2017 ){
      if( !isTightPhotonPOG_Fall17V2( phIdx ) ) return false;
    }
    return true;
    break;

	///////////////
	/// Default ///
	///////////////

  default:
	{
	  cout << "Warning! Photon ID not defined for this id_level!" << endl;
	  return false;
	}

  }//switch
  return true;
}

// Spring 16 ID for 2016 data/MC
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_2016_data_for
// https://indico.cern.ch/event/491548/contributions/2384977/attachments/1377936/2117789/CutBasedPhotonID_25-11-2016.pdf
bool isLoosePhotonPOG_Spring16( int photonIdx )
{
  float pt  = cms3.photons_p4().at(photonIdx).pt();
  float eta = cms3.photons_p4().at(photonIdx).eta();

  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  
  float chiso = photonCHIso03EA(photonIdx, 3);
  float nhiso = photonNHIso03EA(photonIdx, 3);
  float emiso = photonEMIso03EA(photonIdx, 3);

  if(       fabs(eta) < 1.479 ){
    if( hovere > 0.0597                              ) return false;
    if( sieie  > 0.01031                             ) return false;
    if( chiso  > 1.295                               ) return false;
    if( nhiso  > 10.910 + 0.0148*pt + 0.000017*pt*pt ) return false;
    if( emiso  > 3.630 + 0.0047*pt                   ) return false;
  }else if( fabs(eta) > 1.479 ){
    if( hovere > 0.0481                              ) return false;
    if( sieie  > 0.03013                             ) return false;
    if( chiso  > 1.011                               ) return false;
    if( nhiso  > 5.931 + 0.0163*pt + 0.000014*pt*pt  ) return false;
    if( emiso  > 6.641 + 0.0034*pt                   ) return false;
  }
  
  return true;
}

bool isMediumPhotonPOG_Spring16( int photonIdx )
{
  float pt  = cms3.photons_p4().at(photonIdx).pt();
  float eta = cms3.photons_p4().at(photonIdx).eta();

  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  
  float chiso = photonCHIso03EA(photonIdx, 3);
  float nhiso = photonNHIso03EA(photonIdx, 3);
  float emiso = photonEMIso03EA(photonIdx, 3);

  if(       fabs(eta) < 1.479 ){
    if( hovere > 0.0396                              ) return false;
    if( sieie  > 0.01022                             ) return false;
    if( chiso  > 0.441                               ) return false;
    if( nhiso  > 2.725 + 0.0148*pt + 0.000017*pt*pt  ) return false;
    if( emiso  > 2.571 + 0.0047*pt                   ) return false;
  }else if( fabs(eta) > 1.479 ){
    if( hovere > 0.0219                              ) return false;
    if( sieie  > 0.03001                             ) return false;
    if( chiso  > 0.442                               ) return false;
    if( nhiso  > 1.715 + 0.0163*pt + 0.000014*pt*pt  ) return false;
    if( emiso  > 3.863 + 0.0034*pt                   ) return false;
  }
  
  return true;
}

bool isTightPhotonPOG_Spring16( int photonIdx )
{
  float pt  = cms3.photons_p4().at(photonIdx).pt();
  float eta = cms3.photons_p4().at(photonIdx).eta();

  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  
  float chiso = photonCHIso03EA(photonIdx, 3);
  float nhiso = photonNHIso03EA(photonIdx, 3);
  float emiso = photonEMIso03EA(photonIdx, 3);

  if(       fabs(eta) < 1.479 ){
    if( hovere > 0.0269                              ) return false;
    if( sieie  > 0.00994                             ) return false;
    if( chiso  > 0.202                               ) return false;
    if( nhiso  > 0.264 + 0.0148*pt + 0.000017*pt*pt  ) return false;
    if( emiso  > 2.362 + 0.0047*pt                   ) return false;
  }else if( fabs(eta) > 1.479 ){
    if( hovere > 0.0213                              ) return false;
    if( sieie  > 0.03000                             ) return false;
    if( chiso  > 0.034                               ) return false;
    if( nhiso  > 0.586 + 0.0163*pt + 0.000014*pt*pt  ) return false;
    if( emiso  > 2.617 + 0.0034*pt                   ) return false;
  }
  
  return true;
}

// Fall17 V2 ID for 2017+2018 data/MC
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool isLoosePhotonPOG_Fall17V2( int photonIdx )
{
  float pt  = cms3.photons_p4().at(photonIdx).pt();
  float eta = cms3.photons_p4().at(photonIdx).eta();

  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  
  float chiso = photonCHIso03EA(photonIdx, 4);
  float nhiso = photonNHIso03EA(photonIdx, 4);
  float emiso = photonEMIso03EA(photonIdx, 4);

  if(       fabs(eta) < 1.479 ){
    if( hovere > 0.04596                               ) return false;
    if( sieie  > 0.0106                                ) return false;
    if( chiso  > 1.694                                 ) return false;
    if( nhiso  > 24.032 + 0.01512*pt + 2.259e-05*pt*pt ) return false;
    if( emiso  > 2.876 + 0.004017*pt                   ) return false;
  }else if( fabs(eta) > 1.479 ){
    if( hovere > 0.0590                                ) return false;
    if( sieie  > 0.0272                                ) return false;
    if( chiso  > 2.089                                 ) return false;
    if( nhiso  > 19.722 + 0.0117*pt + 2.3e-05*pt*pt    ) return false;
    if( emiso  > 4.162 + 0.0037*pt                     ) return false;
  }
  
  return true;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool isMediumPhotonPOG_Fall17V2( int photonIdx )
{
  float pt  = cms3.photons_p4().at(photonIdx).pt();
  float eta = cms3.photons_p4().at(photonIdx).eta();

  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  
  float chiso = photonCHIso03EA(photonIdx, 4);
  float nhiso = photonNHIso03EA(photonIdx, 4);
  float emiso = photonEMIso03EA(photonIdx, 4);

  if(       fabs(eta) < 1.479 ){
    if( hovere > 0.02197                               ) return false;
    if( sieie  > 0.01015                               ) return false;
    if( chiso  > 1.141                                 ) return false;
    if( nhiso  > 1.189 + 0.01512*pt + 2.259e-05*pt*pt  ) return false;
    if( emiso  > 2.08 + 0.004017*pt                    ) return false;
  }else if( fabs(eta) > 1.479 ){
    if( hovere > 0.0326                                ) return false;
    if( sieie  > 0.0272                                ) return false;
    if( chiso  > 1.051                                 ) return false;
    if( nhiso  > 2.718 + 0.0117*pt + 2.3e-05*pt*pt     ) return false;
    if( emiso  > 3.867 + 0.0037*pt                     ) return false;
  }
  
  return true;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool isTightPhotonPOG_Fall17V2( int photonIdx )
{
  float pt  = cms3.photons_p4().at(photonIdx).pt();
  float eta = cms3.photons_p4().at(photonIdx).eta();

  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  
  float chiso = photonCHIso03EA(photonIdx, 4);
  float nhiso = photonNHIso03EA(photonIdx, 4);
  float emiso = photonEMIso03EA(photonIdx, 4);

  if(       fabs(eta) < 1.479 ){
    if( hovere > 0.02148                               ) return false;
    if( sieie  > 0.00996                               ) return false;
    if( chiso  > 0.65                                  ) return false;
    if( nhiso  > 0.317 + 0.01512*pt + 2.259e-05*pt*pt  ) return false;
    if( emiso  > 2.044 + 0.004017*pt                   ) return false;
  }else if( fabs(eta) > 1.479 ){
    if( hovere > 0.0321                                ) return false;
    if( sieie  > 0.0271                                ) return false;
    if( chiso  > 0.517                                 ) return false;
    if( nhiso  > 2.716 + 0.0117*pt + 2.3e-05*pt*pt     ) return false;
    if( emiso  > 3.032 + 0.0037*pt                     ) return false;
  }
  
  return true;
}

bool isTemplatePhoton( unsigned int phIdx ) {
  if( fabs(photons_p4()            .at(phIdx).pt())  <= 22.0 ) return false; // photon pT > 22 GeV
  if( fabs(photons_p4()            .at(phIdx).eta()) >= 2.4  ) return false; // photon eta < 2.4
  if(photons_full5x5_hOverEtowBC() .at(phIdx)        >= 0.1  ) return false; // hOverE < 0.1
  if(photons_haspixelSeed()        .at(phIdx)                ) return false; // veto pixel seed
  return true;
}

// from https://indico.cern.ch/event/369239/contribution/2/attachments/1134693/1623149/spring15_pcb.pdf
bool isLoosePhoton_Spring15_50ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0103                       ) return false;
	if( chiso  > 2.44                         ) return false;
	if( nhiso  > 2.57 + exp(0.0044*pt+0.5809) ) return false;
	if( emiso  > 1.92 + 0.0043*pt             ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0277                       ) return false;
	if( chiso  > 1.84                         ) return false;
	if( nhiso  > 4.00 + exp(0.0040*pt+0.9402) ) return false;
	if( emiso  > 2.15 + 0.0041*pt             ) return false;
  }
  
  return true;
}

bool isMediumPhoton_Spring15_50ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0100                       ) return false;
	if( chiso  > 1.31                         ) return false;
	if( nhiso  > 0.60 + exp(0.0044*pt+0.5809) ) return false;
	if( emiso  > 1.33 + 0.0043*pt             ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0267                       ) return false;
	if( chiso  > 1.25                         ) return false;
	if( nhiso  > 1.65 + exp(0.0040*pt+0.9402) ) return false;
	if( emiso  > 1.02 + 0.0041*pt             ) return false;
  }
  
  return true;
}

bool isTightPhoton_Spring15_50ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0010                       ) return false;
	if( chiso  > 0.91                         ) return false;
	if( nhiso  > 0.33 + exp(0.0044*pt+0.5809) ) return false;
	if( emiso  > 0.61 + 0.0043*pt             ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0267                       ) return false;
	if( chiso  > 0.65                         ) return false;
	if( nhiso  > 0.93 + exp(0.0040*pt+0.9402) ) return false;
	if( emiso  > 0.54 + 0.0041*pt             ) return false;
  }
  
  return true;
}


// from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Recommended_Working_points_for_2
bool isLoosePhoton_Spring15_25ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                                    ) return false;
	if( sieie  > 0.0102                                  ) return false;
	if( chiso  > 3.32                                    ) return false;
	if( nhiso  > 1.92 + 0.014  * pt + 0.000019*pow(pt,2) ) return false;
	if( emiso  > 0.81 + 0.0053 * pt                      ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere >  0.05                                    ) return false;
	if( sieie  >  0.0274                                  ) return false;
	if( chiso  >  1.97                                    ) return false;
	if( nhiso  > 11.86 + 0.0139 * pt + 0.000025*pow(pt,2) ) return false;
	if( emiso  >  0.83 + 0.0034 * pt                      ) return false;
  }
  
  return true;
}

bool isMediumPhoton_Spring15_25ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                                    ) return false;
	if( sieie  > 0.0102                                  ) return false;
	if( chiso  > 1.37                                    ) return false;
	if( nhiso  > 1.06 + 0.014  * pt + 0.000019*pow(pt,2) ) return false;
	if( emiso  > 0.83 + 0.0034 * pt                      ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                                    ) return false;
	if( sieie  > 0.0268                                  ) return false;
	if( chiso  > 1.10                                    ) return false;
	if( nhiso  > 2.69 + 0.0139 * pt + 0.000025*pow(pt,2) ) return false;
	if( emiso  > 0.39 + 0.0034 * pt                      ) return false;
  }
  
  return true;
}

bool isTightPhoton_Spring15_25ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                                    ) return false;
	if( sieie  > 0.0100                                  ) return false;
	if( chiso  > 0.76                                    ) return false;
	if( nhiso  > 0.97 + 0.014  * pt + 0.000019*pow(pt,2) ) return false;
	if( emiso  > 0.08 + 0.0053 * pt                      ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                                    ) return false;
	if( sieie  > 0.0268                                  ) return false;
	if( chiso  > 0.56                                    ) return false;
	if( nhiso  > 2.09 + 0.0139 * pt + 0.000025*pow(pt,2) ) return false;
	if( emiso  > 0.16 + 0.0034 * pt                      ) return false;
  }
  
  return true;
}

// cuts to emulate trigger cuts from paths: HLT_Photon50_R9Id90_HE10_IsoM etc
bool passTriggerEmu(int photonIdx) {
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();

  if( cms3.photons_full5x5_r9().at(photonIdx) < 0.92     ) return false;
  if( cms3.photons_full5x5_hOverE().at(photonIdx) > 0.12 ) return false;
  if( cms3.photons_tkIsoHollow03().at(photonIdx) > 5     ) return false;
  if( abs(eta) < 1.479 ){
    if( photonEcalpfClusterIso03EA(photonIdx) > 4 + pt * 0.0053  ) return false;
    if( photonHcalpfClusterIso03EA(photonIdx) > 8 + pt * 0.014   ) return false;
  }else if( abs(eta) > 1.479 ){
    if( photonEcalpfClusterIso03EA(photonIdx) > 4 + pt * 0.0034  ) return false;
    if( photonHcalpfClusterIso03EA(photonIdx) > 8 + pt * 0.0139  ) return false;
  }

  return true;
}
