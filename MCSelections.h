#ifndef MCSELECTIONS_H
#define MCSELECTIONS_H

#include "CMS3.h"
#include "TString.h"
#include "Base.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TTree.h"

//Nicely formatted dump of gen-level particles
int  dumpDocLines();

//Parentage functions
bool isFromSUSY(int id, int idx);
bool isFromW(int id, int idx);
bool isFromZ(int id, int idx);
bool isFromB(int id, int idx);
bool isFromC(int id, int idx);
bool isFromLight(int id, int idx);
bool isFromLightFake(int id, int idx);
bool idIsCharm(int id);
bool idIsBeauty(int id);

// sourceId to match CERN/ETH
int getSourceId(int genpIdx);
bool useSourceId(int sourceId);
bool terminateSourceId(int motherId);
bool terminateMotherId(int motherId);

// top pt reweighting
float topPtWeight(float pt_top, float pt_tbar);
float topPtWeight_cutoff(float pt_top, float pt_tbar);

// 2016 ISR recipe
int get_nisrMatch ( const std::vector<LorentzVector>& clean_jets );
float get_isrWeight ( int nisrMatch );
float get_isrUnc ( int nisrMatch );
float get_isrWeight ( int nisrMatch, int year );
float get_isrUnc ( int nisrMatch, int year );
// Moriond 2017 ISR recipe for ewkino models
float get_isrWeight_ewk ( float isr_pt );

#endif
