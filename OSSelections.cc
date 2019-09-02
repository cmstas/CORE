#include <stdexcept>

#include "OSSelections.h"
#include "MuonSelections.h"
#include "ElectronSelections.h"
#include "PhotonSelections.h"
#include "JetSelections.h"
#include "VertexSelections.h"
#include "Math/VectorUtil.h"
#include "Tools/utils.h"

using namespace tas;

float els_dzPV_firstPV(unsigned int elIdx)
{
/*This function exists because nanoAOD does not have any quality checks
for PV, and computes dZ and dxy blindly with the first PV*/

    float dz = -999;

    int first_good_vertex = firstGoodVertex();
    if(first_good_vertex == 0)
        return els_dzPV().at(elIdx);

    LorentzVector track_vertex = els_trk_vertex_p4().at(elIdx);
    LorentzVector primary_vertex = vtxs_position().at(0);
    LorentzVector track_p4 = els_trk_p4().at(elIdx);

        //Reproducing dz computation from CMSSW : https://github.com/cms-sw/cmssw/blob/0b6f294b44aef89ab2b5251053531bbdb212f720/DataFormats/TrackReco/interface/TrackBase.h#L666
    dz = (track_vertex.Z() - primary_vertex.Z()) - ((track_vertex.X() - primary_vertex.X()) * track_p4.X() + (track_vertex.Y() - primary_vertex.Y()) * track_p4.Y())/track_p4.Pt() * track_p4.Z()/track_p4.Pt();

    return dz;    
}


float els_dxyPV_firstPV(unsigned int elIdx)
{
    /*This function exists because nanoAOD does not have any quality checks
for PV, and computes dZ and dxy blindly with the first PV*/ 

    float dxy = -999;
    int first_good_vertex = firstGoodVertex();
    if(first_good_vertex == 0)
        return els_dxyPV().at(elIdx);

    LorentzVector track_vertex = els_trk_vertex_p4().at(elIdx);
    LorentzVector primary_vertex = vtxs_position().at(0);
    LorentzVector track_p4 = els_trk_p4().at(elIdx);

    //Reproducing dz computation from CMSSW : https://github.com/cms-sw/cmssw/blob/0b6f294b44aef89ab2b5251053531bbdb212f720/DataFormats/TrackReco/interface/TrackBase.h#L646
   
    dxy =  (-(track_vertex.X() - primary_vertex.X()) * track_p4.Y() + (track_vertex.Y() - primary_vertex.Y()) * track_p4.X()) / track_p4.Pt();

        return dxy;
}

//Since we don't have muons track vertices stored in CMS4, we go to the isotracks collection to get them muon tracks out
int mus_findOverlapIsotrack(unsigned int muIdx)
{
   float lowestDR = 0.1;
   size_t closestIsotrack = -1;
   for(size_t iit = 0; iit < cms3.isotracks_p4().size(); iit++)
   {
       float currentDR = DeltaR(cms3.isotracks_p4().at(iit),cms3.mus_p4().at(muIdx));  
	if(currentDR < lowestDR && abs(cms3.isotracks_particleId().at(iit)) == 13)
       {
	    lowestDR = currentDR;
            closestIsotrack = iit;
       }
   }
   return closestIsotrack;
}

float mus_dzPV_firstPV(unsigned int muIdx)
{
    int first_good_vertex = firstGoodVertex();
    float dz_difference = 0;
    if(first_good_vertex == 0)
        return mus_dzPV().at(muIdx);

    //Fancy ass dz calculation comin' up
   
    LorentzVector best_PV = vtxs_position().at(first_good_vertex);
    LorentzVector first_PV = vtxs_position().at(0);
    LorentzVector track_p4 = mus_trk_p4().at(muIdx);

    //Stop-gap for 0 inner track Pt
    if(track_p4.pt() == 0)
        return -9999; 

    dz_difference = (best_PV.Z() - first_PV.Z()) - ((best_PV.X() - first_PV.X()) * track_p4.X() + (best_PV.Y() - first_PV.Y()) * track_p4.Y())/track_p4.Pt() * track_p4.Z()/track_p4.Pt();

    return mus_dzPV().at(muIdx) + dz_difference;    
}

/*float mus_dzPV_firstPV(unsigned int muIdx)
{
    int first_good_vertex = firstGoodVertex();
    if(first_good_vertex == 0)
        return mus_dzPV().at(muIdx);
    
    int overlap_isotrack_index = mus_findOverlapIsotrack(muIdx);
    if(overlap_isotrack_index >= 0)

    {
        return cms3.isotracks_dz().at(overlap_isotrack_index); 
    }
    else
        return -999;
}*/

float mus_dxyPV_firstPV(unsigned int muIdx)
{
    int first_good_vertex = firstGoodVertex();
    float dxy_difference = 0;
    if(first_good_vertex == 0)
        return mus_dxyPV().at(muIdx);

    LorentzVector best_PV = vtxs_position().at(first_good_vertex);
    LorentzVector first_PV = vtxs_position().at(0);
    LorentzVector track_p4 = mus_trk_p4().at(muIdx);

    if(track_p4.pt() == 0)
        return -9999;

    dxy_difference = (-(best_PV.X() - first_PV.X()) * track_p4.Y() + (best_PV.Y() - first_PV.Y()) * track_p4.X()) / track_p4.Pt();

    return mus_dxyPV().at(muIdx) + dxy_difference;


}

/*float mus_dxyPV_firstPV(unsigned int muIdx)
{
    int first_good_vertex = firstGoodVertex();
    if(first_good_vertex == 0)
        return mus_dxyPV().at(muIdx);

    int overlap_isotrack_index = mus_findOverlapIsotrack(muIdx);
    if(overlap_isotrack_index >= 0)
    {
	    return cms3.isotracks_dxy().at(overlap_isotrack_index); 
    }
    else
        return -999;

}*/


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
  return false;
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
  if (abs(els_dzPV_firstPV(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV_firstPV(index)) >= 0.05                      ) return false;// dR < 0.05
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
  if (abs(els_dzPV_firstPV(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV_firstPV(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  return true;
}

bool passElectronSelection_ZMET_v8(int index, bool vetoTransition, bool eta24)
{
    //2018 stuff when it comes out
    if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
    if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
    if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
    if( !electronID( index, ZMET_tightMVA_v4 )       ) return false; // Electron ID  

    //IP & trigger cuts to be compatible with multilepton baseline cuts
    if (abs(els_dzPV_firstPV(index)) >= 0.1                       ) return false;// dZ < 0.1
    if (abs(els_dxyPV_firstPV(index)) >= 0.05                      ) return false;// dR < 0.05
    if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
    return true;

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
    if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 10 GeV - Minimum pT cut
    if( vetoTransition
	  && fabs(cms3.els_p4().at(index).eta()) > 1.4
	  && fabs(cms3.els_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
    if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.5    ) return false; // eta < 2.5


      //IP & trigger cuts to be compatible with multilepton baseline cuts
    if (abs(els_dzPV_firstPV(index)) >= 0.1                       ) return false;// dZ < 0.1
    if (abs(els_dxyPV_firstPV(index)) >= 0.05                      ) return false;// dR < 0.05
    if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8

    if( !electronID( index, ZMET_looseMVA_v3 )       ) return false; // Electron ID  

    return true;
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
  if (abs(els_dzPV_firstPV(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV_firstPV(index)) >= 0.05                      ) return false;// dR < 0.05
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
  if (abs(els_dzPV_firstPV(index)) >= 0.1                       ) return false;// dZ < 0.1
  if (abs(els_dxyPV_firstPV(index)) >= 0.05                      ) return false;// dR < 0.05
  if (abs(els_ip3d()  .at(index))/els_ip3derr().at(index) >= 8 ) return false;// SIP3D < 8
  return true;
}


//~-~-~-~-~-~-~-~//
//Muon selections//
//~-~-~-~-~-~-~-~//
bool passMuonSelection_ZMET(int index){

  int year = gconf.year;
  if(year == 2016)
    return passMuonSelection_ZMET_v7( index, true, true );
  else if(year == 2017)
      return passMuonSelection_ZMET_v8(index,true,true);
  else if(year==2018)
      return passMuonSelection_ZMET_v9(index,true,true);
  return false;
}


bool passMuonSelection_ZMET_v9(int index, bool vetoTransition, bool eta24)
{
    //2018 data - update when it comes
    if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
    if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
    if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
    if( !muonID( index, ZMET_mediumMu_v4 )              ) return false; // medium Muon ID  

    //IP cuts to be compatible with multilepton baseline cuts
    if (fabs(mus_dxyPV_firstPV(index))   > 0.05 ) return false;
	if (fabs(mus_dzPV_firstPV(index))   > 0.1  ) return false;
    if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
    return true;

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
  if (fabs(mus_dxyPV_firstPV(index))   > 0.05 ) return false;
  if (fabs(mus_dzPV_firstPV(index))   > 0.1  ) return false;
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
  if (fabs(mus_dxyPV_firstPV(index))   > 0.05 ) return false;
  if (fabs(mus_dzPV_firstPV(index))   > 0.1  ) return false;
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
  if (fabs(mus_dxyPV_firstPV(index))   > 0.05 ) return false;
  if (fabs(mus_dzPV_firstPV(index))   > 0.1  ) return false;
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


bool passMuonSelection_ZMET_veto(int index)
{
    if(gconf.year == 2016)
    {
        passMuonSelection_ZMET_veto_v3(index,false,true);
    }
    else if(gconf.year == 2017 || gconf.year == 2018)
    {
        passMuonSelection_ZMET_veto_v4(index,false,true);
    }

    return false;
}


bool passMuonSelection_ZMET_veto_v4(int index, bool vetoTransition, bool eta24)
{
    if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.6  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_veto_v4 )         ) return false; // medium Muon ID with looser iso

  //IP cuts to be compatible with multilepton baseline cuts
  if (fabs(mus_dxyPV() .at(index))   > 0.05 ) return false;
  if (fabs(mus_dzPV()  .at(index))   > 0.1  ) return false;
  if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
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
  if (fabs(mus_dxyPV() .at(index))   > 0.05 ) return false;
  if (fabs(mus_dzPV()  .at(index))   > 0.1  ) return false;
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
  if (fabs(mus_dxyPV() .at(index))   > 0.05 ) return false;
  if (fabs(mus_dzPV()  .at(index))   > 0.1  ) return false;

  if (abs(mus_ip3d().at(index))/mus_ip3derr().at(index) >= 8) return false;// sip3d < 8
  return true;
}

bool passPhotonSelection_ZMET(int index ){
    if(gconf.year == 2016)
        return passPhotonSelection_ZMET_v4(index, true, true);
    else if(gconf.year == 2017)
        return passPhotonSelection_ZMET_v5(index, true, true);
    else if(gconf.year == 2018)
        return passPhotonSelection_ZMET_v5(index, true, true); //2017 and 2018 same ID 
    return false;
}

bool passPhotonSelection_ZMET_v6(int index, bool vetoTransition, bool eta24)
{
    //2018 not yet here
    return false;
}

bool passPhotonSelection_ZMET_v5(int index, bool vetoTransition, bool eta24)
{
  if( fabs(cms3.photons_p4().at(index).pt()) < 22.0       ) return false; // pT > 50 GeV - Minimum pT cut based on 2017 triggers
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapElectron_ZMET_v2( index, 10.0 )              ) return false; // remove electrons from W
  if(!photonID(index,ZMET_photon_v5)) return false;
  return true;
 
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

