#include "CMS3.h"
#include "Base.h"

#ifndef OSSELECTIONS_H
#define OSSELECTIONS_H


float els_dzPV_firstPV(unsigned int elIdx);
float els_dxyPV_firstPV(unsigned int elIdx);
float mus_dzPV_firstPV(unsigned int muIdx);
float mus_dxyPV_firstPV(unsigned int muIdx);
int mus_findOverlapIsotrack(unsigned int muIdx);

bool overlapMuon_ZMET_v1    ( int index , float ptcut = 10.0 );
bool overlapElectron_ZMET_v2( int index , float ptcut = 10.0);
bool overlapElectron_ZMET_v1( int index , float ptcut = 10.0);

bool passElectronSelection_ZMET(int index);
bool passElectronSelection_ZMET_veto(int index );
bool passElectronSelection_ZMET_v8(int index, bool vetoTransition, bool eta24);
bool passElectronSelection_ZMET_v7(int index, bool vetoTransition, bool eta24);
bool passElectronSelection_ZMET_v6(            int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_v5(            int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_v4(            int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_v3(            int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_v2(            int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_NoIso_v2(      int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_v1(            int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_v1_NoIso(      int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_thirdlepton_v1(int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_thirdlepton_v2(int index, bool vetoTransition, bool eta24 );
bool passElectronSelection_ZMET_thirdlepton_v3(int index, bool vetoTransition, bool eta24);
bool passElectronSelection_ZMET_thirdlepton_v4(int index, bool vetoTransition, bool eta24);


bool passMuonSelection_ZMET(int index);
bool passMuonSelection_ZMET_v9(int index, bool vetoTransition, bool eta24);
bool passMuonSelection_ZMET_v8(int index, bool vetoTransition, bool eta24);
bool passMuonSelection_ZMET_v7(       int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_v6(       int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_v5(       int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_v4(       int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_v3(       int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_v2(       int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_NoIso_v2( int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_v1_NoIso( int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_v1(       int index, bool vetoTransition, bool eta24 );

bool passMuonSelection_ZMET_veto(int index);
bool passMuonSelection_ZMET_veto_v4(int index, bool vetoTransition, bool eta24);
bool passMuonSelection_ZMET_veto_v3(  int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_veto_v2(  int index, bool vetoTransition, bool eta24 );
bool passMuonSelection_ZMET_veto_v1(  int index, bool vetoTransition, bool eta24 );

bool passPhotonSelection_ZMET(int index );
bool passPhotonSelection_ZMET_v6(int index, bool vetoTransition, bool eta24 );
bool passPhotonSelection_ZMET_v5(int index, bool vetoTransition, bool eta24 );
bool passPhotonSelection_ZMET_v4(int index, bool vetoTransition, bool eta24 );
bool passPhotonSelection_ZMET_v3(int index, bool vetoTransition, bool eta24 );
bool passPhotonSelection_ZMET_v2(int index, bool vetoTransition, bool eta24 );
bool passPhotonSelection_ZMET_v1(int index, bool vetoTransition, bool eta24 );

bool electronPassesHLTEmulator(int index);

float mbb_highest_csv(std::vector <LorentzVector> jets_p4, std::vector<float> jets_csv);

bool passJetSelection_ZMET(int iJet);
bool passBTagWP_ZMET(int iJet, int bTagIndex, int WP);

#endif
