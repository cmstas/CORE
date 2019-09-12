#include <stdexcept>

#include "VVVSelections.h"
#include "MuonSelections.h"
#include "ElectronSelections.h"
#include "IsolationTools.h"

#include "Math/VectorUtil.h"

using namespace tas;

/* ID Levels:
VVV_cutbased_veto
VVV_cutbased_veto_noiso
VVV_cutbased_veto_noiso_noip
VVV_cutbased_veto_bak1             //electrons only
VVV_cutbased_veto_noiso_bak1       //electrons only
VVV_cutbased_veto_noiso_noip_bak1  //electrons only
VVV_cutbased_fo
VVV_cutbased_fo_noiso
VVV_cutbased_tight_noiso
VVV_cutbased_tight
VVV_MVAbased_tight_noiso           //electrons only
VVV_MVAbased_tight                 //electrons only
VVV_baseline
*/

//~-~-~-~-~-~-~-~-~-~//
//Electron selections//
//~-~-~-~-~-~-~-~-~-~//
bool passElectronSelection_VVV(int index, id_level_t id_string){
  return passElectronSelection_VVV_v1( index, id_string, true, true );
}

bool passElectronSelection_VVV_v1(int index, id_level_t id_string, bool vetoTransition, bool eta24){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0 ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition && fabs(cms3.els_p4().at(index).eta()) > 1.4 && fabs(cms3.els_p4().at(index).eta()) < 1.6 ) return false; // veto x-ition region
  if( eta24 && fabs(cms3.els_p4()[index].eta()) > 2.4 ) return false; // eta < 2.5
  if( !electronID( index, id_string ) ) return false; // Electron ID  
  return true;
}

//~-~-~-~-~-~-~-~//
//Muon selections//
//~-~-~-~-~-~-~-~//
bool passMuonSelection_VVV(int index, id_level_t id_string){
  return passMuonSelection_VVV_v1( index, id_string, false, true );
}

bool passMuonSelection_VVV_v1(int index, id_level_t id_string, bool vetoTransition, bool eta24){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0 ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition && fabs(cms3.mus_p4().at(index).eta()) > 1.4 && fabs(cms3.mus_p4().at(index).eta()) < 1.6 ) return false; // veto x-ition region
  if( eta24 && fabs(cms3.mus_p4()[index].eta()) > 2.4 ) return false; // eta < 2.5
  if( !muonID( index, id_string ) ) return false; // Muon ID  
  return true;
}

//##############################################################################################################
// Very Loose Lepton ID
bool isPt10POGVetoElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 10.          )) return false;
    if (!( isVetoElectronPOGfall17_v2(idx)        )) return false;
    if (!( fabs(cms3.els_p4()[idx].eta()) < 2.5   )) return false;
    if (fabs(cms3.els_etaSC()[idx]) <= 1.479)
    {
        if (!( fabs(cms3.els_dzPV()[idx]) < 0.1       )) return false;
        if (!( fabs(cms3.els_dxyPV()[idx]) < 0.05     )) return false;
    }
    else
    {
        if (!( fabs(cms3.els_dzPV()[idx]) < 0.2       )) return false;
        if (!( fabs(cms3.els_dxyPV()[idx]) < 0.1      )) return false;
    }
    if (!( fabs(cms3.els_ip3d()[idx] / cms3.els_ip3derr()[idx]) < 4. )) return false;
    return true;
}

//##############################################################################################################
// Very Loose Lepton ID
bool isPt10POGVetoMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 10.        )) return false;
    if (!( isMediumMuonPOG(idx)                 )) return false;
    if (!( fabs(cms3.mus_p4()[idx].eta()) < 2.4 )) return false;
    if (!( muRelIso04DB(idx)  < 0.25            )) return false;
    if (fabs(cms3.mus_p4()[idx].eta()) <= 1.479)
    {
        if (!( fabs(cms3.mus_dzPV()[idx]) < 0.1       )) return false;
        if (!( fabs(cms3.mus_dxyPV()[idx]) < 0.05     )) return false;
    }
    else
    {
        if (!( fabs(cms3.mus_dzPV()[idx]) < 0.2       )) return false;
        if (!( fabs(cms3.mus_dxyPV()[idx]) < 0.1      )) return false;
    }
    if (!( fabs(cms3.mus_ip3d()[idx] / cms3.mus_ip3derr()[idx]) < 4. )) return false;
    return true;
}
