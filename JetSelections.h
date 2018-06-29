#ifndef JETSELECTIONS_H
#define JETSELECTIONS_H
#include "CMS3.h"
#include "TString.h"
#include "Base.h"

// recommended jet id for 2017 94x analyses
bool isTightPFJet_2017_v1(unsigned int pfJetIdx);

// Jet ID as of Sept. 11th based on this twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetID?rev=95
bool isLoosePFJet_Summer16_v1(unsigned int pfJetIdx, bool use_puppi = false);
bool isTightPFJet_Summer16_v1(unsigned int pfJetIdx);
bool isTightPFJetLepVeto_Summer16_v1(unsigned int pfJetIdx);

bool isLoosePFJet(unsigned int pfJetIdx);
bool isMediumPFJet(unsigned int pfJetIdx);
bool isTightPFJet(unsigned int pfJetIdx);

bool isTightPFJet2017(unsigned int pfJetIdx);

bool isLoosePFJetV2(unsigned int pfJetIdx);
bool isTightPFJetV2(unsigned int pfJetIdx);

bool isLoosePFJet_50nsV1(unsigned int pfJetIdx, bool use_puppi = false);
bool isTightPFJet_50nsV1(unsigned int pfJetIdx);
bool isTightPFJetLepVeto_50nsV1(unsigned int pfJetIdx);
bool isMonoPFJet_MT2(unsigned int pfJetIdx);
bool isMonoPFJet_Monojet(unsigned int pfJetIdx);

bool loosePileupJetId(unsigned int pfJetIdx);
bool loosePileupJetId_v2(unsigned int pfJetIdx, bool use_puppi = false);

bool JetIsElectron(LorentzVector pfJet, id_level_t id_level, float ptcut = 7., float deltaR = 0.4);
bool JetIsMuon(LorentzVector pfJet, id_level_t id_level, float ptcut = 5., float deltaR = 0.4);

bool isBadFastsimJet(unsigned int pfJetIdx);

#endif
