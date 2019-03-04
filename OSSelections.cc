#include <stdexcept>

#include "OSSelections.h"
#include "MuonSelections.h"
#include "ElectronSelections.h"
#include "PhotonSelections.h"
#include "JetSelections.h"

#include "Math/VectorUtil.h"

using namespace tas;

bool overlapMuon_ZMET_v1( int index , float ptcut ){
  for( unsigned int muind = 0; muind < cms3.mus_p4().size(); muind++ ){
    float dr = ROOT::Math::VectorUtil::DeltaR( cms3.els_p4().at(index) , cms3.mus_p4().at(muind) );    
    if( dr > 0.05                                     ) continue;
    if( cms3.mus_p4().at(muind).pt() < ptcut          ) continue;
	if( !passMuonSelection_ZMET_v2(muind, true, true) ) continue;
    return true;
  }
  return false;
}

bool overlapElectron_ZMET_v3(int index, float ptcut)
{
    for( unsigned int elind = 0; elind < cms3.els_p4().size(); elind++ ){
    float dr = ROOT::Math::VectorUtil::DeltaR( cms3.photons_p4().at(index) , cms3.els_p4().at(elind) );    
    if( dr > 0.2                                          ) continue;
    if( cms3.els_p4().at(elind).pt() < ptcut              ) continue;
    if( !passElectronSelection_ZMET_veto(elind) ) continue;
    return true;
  }
  return false;

}

bool overlapElectron_ZMET_v2( int index , float ptcut ){
  for( unsigned int elind = 0; elind < cms3.els_p4().size(); elind++ ){
    float dr = ROOT::Math::VectorUtil::DeltaR( cms3.photons_p4().at(index) , cms3.els_p4().at(elind) );    
    if( dr > 0.2                                          ) continue;
    if( cms3.els_p4().at(elind).pt() < ptcut              ) continue;
    if( !passElectronSelection_ZMET_veto(elind) ) continue;
    return true;
  }
  return false;
}

bool overlapElectron_ZMET_v1( int index , float ptcut ){
  for( unsigned int elind = 0; elind < cms3.els_p4().size(); elind++ ){
    float dr = ROOT::Math::VectorUtil::DeltaR( cms3.photons_p4().at(index) , cms3.els_p4().at(elind) );    
    if( dr > 0.2                                          ) continue;
    if( cms3.els_p4().at(elind).pt() < ptcut              ) continue;
	if( !passElectronSelection_ZMET_v2(elind, true, true) ) continue;
    return true;
  }
  return false;
}

//~-~-~-~-~-~-~-~-~-~//
//Electron selections//
//~-~-~-~-~-~-~-~-~-~//
bool passElectronSelection_ZMET(int index){

  int year = gconf.year;
  if(year == 2016)
    return passElectronSelection_ZMET_v6( index, true, true );
  else if(year == 2017)
      return passElectronSelection_ZMET_v7(index, true, true);
  else if(year == 2018)
      return passElectronSelection_ZMET_v8(index, true, true);
  return false;
}

bool passElectronSelection_ZMET_veto(int index ){
  if(gconf.year == 2016)
    return passElectronSelection_ZMET_thirdlepton_v2( index, false, false );
  else if(gconf.year == 2017)
      return passElectronSelection_ZMET_thirdlepton_v3(index, false, false);
  else if(gconf.year == 2018)
      return passElectronSelection_ZMET_thirdlepton_v4(index, false, false);
}

bool passElectronSelection_ZMET_v6(int index, bool vetoTransition, bool eta24 ){
    //2016 stuff
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  if( !electronID( index, ZMET_tightMVA_v2 )       ) return false; // Electron ID  

  //IP & trigger cuts to be compatible with multilepton baseline cuts
  if (abs(els_dzPV()  .at(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV() .at(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  return true;
}


bool passElectronSelection_ZMET_v7(int index, bool vetoTransition, bool eta24)
{
    //2017 stuff
    if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  if( !electronID( index, ZMET_tightMVA_v3 )       ) return false; // Electron ID  

  //IP & trigger cuts to be compatible with multilepton baseline cuts
  if (abs(els_dzPV()  .at(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV() .at(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  return true;
}

bool passElectronSelection_ZMET_v8(int index, bool vetoTransition, bool eta24)
{
    //2018 stuff when it comes out
    return false;
}

bool passElectronSelection_ZMET_v5(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  if( !electronID( index, ZMET_tightMVA_v2 )       ) return false; // Electron ID  

  //IP & trigger cuts to be compatible with multilepton baseline cuts
  if (abs(els_dzPV()  .at(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV() .at(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  if( !electronPassesHLTEmulator(index)                        ) return false;// emulate trigger cuts
  return true;
}

bool passElectronSelection_ZMET_v4(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_tightMVA_v2 )       ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_v3(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_tightMVA_v1 )       ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_v2 )          ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_NoIso_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_noiso_v2 )    ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_v1_NoIso(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMet2012_v1(index,10.0)          ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_noiso_v1 )    ) return false; // Electron ID  
  return true;
}


bool passElectronSelection_ZMET_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMet2012_v1(index,10.0)          ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_v1 )          ) return false; // Electron ID  
  return true;
}

//This is the 3rd lepton veto ID to be in sync with the edge/multilepton group
bool passElectronSelection_ZMET_thirdlepton_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.5    ) return false; // eta < 2.5
  if( !electronID( index, ZMET_looseMVA_v1 )       ) return false; // Electron ID  

  //IP & trigger cuts to be compatible with multilepton baseline cuts
  if (abs(els_dzPV()  .at(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV() .at(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  if( !electronPassesHLTEmulator(index)                        ) return false;// emulate trigger cuts
  return true;
}


bool passElectronSelection_ZMET_thirdlepton_v4(int index, bool vetoTransition, bool eta24)
{
    return false;
}

bool passElectronSelection_ZMET_thirdlepton_v3(int index, bool vetoTransition, bool eta24)
{
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.5    ) return false; // eta < 2.5
  if( !electronID( index, ZMET_looseMVA_v2 )       ) return false; // Electron ID  

  //IP & trigger cuts to be compatible with multilepton baseline cuts
  if (abs(els_dzPV()  .at(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV() .at(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  return true;

}

//This is the 3rd lepton veto ID to be in sync with the edge/multilepton group, trigger emulation removed
bool passElectronSelection_ZMET_thirdlepton_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.5    ) return false; // eta < 2.5
  if( !electronID( index, ZMET_looseMVA_v1 )       ) return false; // Electron ID  

  //IP & trigger cuts to be compatible with multilepton baseline cuts
  if (abs(els_dzPV()  .at(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV() .at(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  return true;
}


//~-~-~-~-~-~-~-~//
//Muon selections//
//~-~-~-~-~-~-~-~//
bool passMuonSelection_ZMET(int index){

  int year = gconf.year;
  if(year == 2016)
    return passMuonSelection_ZMET_v5( index, true, true );
  else if(year == 2017)
      return passMuonSelection_ZMET_v8(index,true,true);
  else if(year==2018)
      return passMuonSelection_ZMET_v9(index,true,true);
  return false;
}


bool passMuonSelection_ZMET_v9(int index, bool vetoTransition, bool eta24)
{
    //2018 data - update when it comes
    return false;
}

bool passMuonSelection_ZMET_v8(int index, bool vetoTransition, bool eta24)
{
  //2017 data and MC

  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_v4 )              ) return false; // medium Muon ID  

  //IP cuts to be compatible with multilepton baseline cuts
  if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
  return true;
 
}


//278820 start or 2016G
bool passMuonSelection_ZMET_v7(int index, bool vetoTransition, bool eta24 ){
  if( cms3.evt_isRealData() ){
	if( cms3.evt_run() >= 278820 ) return passMuonSelection_ZMET_v5( index, true, true );
	if( cms3.evt_run() <  278820 ) return passMuonSelection_ZMET_v6( index, true, true );
  }
  else{
	return passMuonSelection_ZMET_v5( index, true, true );
  }
  cout<<"Warning! Muon ID selection should not get here. Please check the selection."<<endl;
  return false;
}

bool passMuonSelection_ZMET_v6(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_v3 )              ) return false; // medium Muon ID  

  //IP cuts to be compatible with multilepton baseline cuts
  if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
  return true;
}

bool passMuonSelection_ZMET_v5(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_v2 )              ) return false; // medium Muon ID  

  //IP cuts to be compatible with multilepton baseline cuts
  if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
  return true;
}

bool passMuonSelection_ZMET_v4(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_v2 )              ) return false; // medium Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v3(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 15.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
  	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
  	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_v1 )              ) return false; // medium Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_v2 )                 ) return false; // tight Muon ID  
  return true;
}

bool passMuonSelection_ZMET_NoIso_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 15.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_noiso_v2 )           ) return false; // tight Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v1_NoIso(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_noiso_v1 )           ) return false; // tight Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_v1 )                 ) return false; // tight Muon ID  
  return true;
}

//278820 start or 2016G
bool passMuonSelection_ZMET_veto_v3(int index, bool vetoTransition, bool eta24 ){
  if( cms3.evt_isRealData() ){
	if( cms3.evt_run() >= 278820 ) return passMuonSelection_ZMET_veto_v1( index, false, true );
	if( cms3.evt_run() <  278820 ) return passMuonSelection_ZMET_veto_v2( index, false, true );
  }
  else{
	return passMuonSelection_ZMET_veto_v1( index, false, true );
  }
  cout<<"Warning! Muon ID selection should not get here. Please check the selection."<<endl;
  return false;
}

  // veto selection to be used on 3rd lepton to compare with edge/multilepton group
bool passMuonSelection_ZMET_veto_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_veto_v3 )         ) return false; // medium Muon ID with looser iso

  //IP cuts to be compatible with multilepton baseline cuts
  if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
  return true;
}

// veto selection to be used on 3rd lepton to compare with edge/multilepton group
bool passMuonSelection_ZMET_veto_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_veto_v2 )         ) return false; // medium Muon ID with looser iso

  //IP cuts to be compatible with multilepton baseline cuts
  if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
  return true;
}

bool passPhotonSelection_ZMET(int index ){
    if(gconf.year == 2016)
        return passPhotonSelection_ZMET_v4(index, true, true);
    else if(gconf.year == 2017)
        return passPhotonSelection_ZMET_v5(index, true, true);
    else if(gconf.year == 2018)
        return passPhotonSelection_ZMET_v6(index, true, true);
}

bool passPhotonSelection_ZMET_v6(int index, bool vetoTransition, bool eta24)
{
    //2018 not yet here
    return false;
}

bool passPhotonSelection_ZMET_v5(int index, bool vetoTransition, bool eta24)
{
   if( fabs(cms3.photons_p4().at(index).pt()) < 50.0       ) return false; // pT > 50 GeV - Minimum pT cut based on 2017 triggers
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapElectron_ZMET_v2( index, 10.0 )              ) return false; // remove electrons from W
 
}

// move trigger emulation cuts on pfcluster iso etc into PhotonSelections in CORE
bool passPhotonSelection_ZMET_v4(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.photons_p4().at(index).pt()) < 50.0       ) return false; // pT > 50 GeV - Minimum pT cut based on 2017 triggers
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapElectron_ZMET_v2( index, 10.0 )              ) return false; // remove electrons from W
  if( !photonID(index, ZMET_photon_v4 )                   ) return false;
  return true;
}

  bool passPhotonSelection_ZMET_v3(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.photons_p4().at(index).pt()) < 22.0       ) return false; // pT > 22 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapElectron_ZMET_v1( index, 10.0 )              ) return false; // remove electrons from W
  if( !photonID(index, ZMET_photon_v3 )                   ) return false;
  return true;
}

  bool passPhotonSelection_ZMET_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.photons_p4().at(index).pt()) < 22.0       ) return false; // pT > 22 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapElectron_ZMET_v1( index, 10.0 )              ) return false; // remove electrons from W
  if( !photonID(index, ZMET_photon_v1 )                   ) return false;
  return true;
}

bool passPhotonSelection_ZMET_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.photons_p4().at(index).pt()) < 22.0       ) return false; // pT > 22 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !photonID(index, ZMET_photon_v1 ) ) return false;
  return true;
}

// used by the multilepton group since there are no triggers in 80X
bool electronPassesHLTEmulator(int index){

  float etasc = cms3.els_etaSC().at(index);
  if (    cms3.els_hOverE().at(index)  >= (0.10 - 0.03  * (abs(etasc)>1.479))) return false;
  if (abs(cms3.els_dEtaIn().at(index)) >= (0.01 - 0.002 * (abs(etasc)>1.479))) return false;
  if (abs(cms3.els_dPhiIn().at(index)) >= (0.04 + 0.03  * (abs(etasc)>1.479))) return false;

  float eInvMinusPInv = (1.0/els_ecalEnergy().at(index)) - (els_eOverPIn().at(index)/els_ecalEnergy().at(index));
  if (eInvMinusPInv <= -0.05                           ) return false;
  if (eInvMinusPInv >= (0.01-0.005*(abs(etasc)>1.479)) ) return false;

  if (cms3.els_sigmaIEtaIEta_full5x5().at(index) >= (0.011+0.019*(abs(etasc)>1.479)) ) return false;

  return true;
}

// This is meant to be passed as the third argument, the predicate, of the standard library sort algorithm
inline bool sortbyCSV( pair<LorentzVector, float> &vec1, pair<LorentzVector, float> &vec2 ) {
  return vec1.second > vec2.second;
}


float mbb_highest_csv(vector <LorentzVector> jets_p4, vector<float> jets_csv){

  float mbb = 0.0;

  if( jets_p4.size() < 2 ){ mbb = -99.0;}
  else{

	vector <pair<LorentzVector, float>> jet_csv_pairs;
	for( size_t jetind = 0; jetind < jets_p4.size(); jetind++ ){
	  if( jets_csv.at(jetind) > 0.8484){
		jet_csv_pairs.push_back(make_pair(jets_p4.at(jetind),jets_csv.at(jetind)));
	  }
	}

	sort( jet_csv_pairs.begin(), jet_csv_pairs.end(), sortbyCSV);

	if( jet_csv_pairs.size() > 1 ){
	  mbb = (jet_csv_pairs.at(0).first + jet_csv_pairs.at(1).first).mass();
	}else{
	  mbb = -99.0;
	}
  
  }
  return mbb;
}


//JET SELECTIONS
bool passJetSelection_ZMET(int iJet)
{
    if(gconf.year == 2016)
        return isLoosePFJet_Summer16_v1(iJet);
    else if(gconf.year == 2017)
        return isTightPFJet_2017_v1(iJet);
    else if(gconf.year == 2018)
        return false; //Not out yet!
    return false;
}

//B-TAG SELECTIONS
bool passBTagWP_ZMET(int iJet, int bTagIndex, int WP)
{
    /*
    * bTagIndex - 1 for DeepCSV, more will be added later if required
    * WP - 0 for loose, 1 for medium, 2 for tight
    */    

    float looseWP = 0.0, mediumWP = 0.0, tightWP = 0.0;
    
    if(gconf.year == 2016)
    {
        if(bTagIndex == 1)
        {
               
        }
        else if(bTagIndex == 2)
        {

        }
    }
    else if(gconf.year == 2017)
    {
        if(bTagIndex == 1)
        {
            looseWP = 0.1522;
            mediumWP = 0.4941;
            tightWP = 0.8001;
        }
    }
    else if(gconf.year == 2018)
    {
        //not released yet
    }


    return false;
    
}

