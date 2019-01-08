#ifndef PHOTONSELECTIONS_H
#define PHOTONSELECTIONS_H
#include "CMS3.h"
#include "TString.h"
#include "Base.h"

bool isLoosePhoton(unsigned int phIdx, analysis_t analysis, int version = 1);
bool isTightPhoton(unsigned int phIdx, analysis_t analysis, int version = 1);
bool photonID(unsigned int phIdx, id_level_t id_level);
bool isTemplatePhoton( unsigned int phIdx );
bool isLoosePhoton_Spring15_50ns( int photonIdx );
bool isMediumPhoton_Spring15_50ns( int photonIdx );
bool isTightPhoton_Spring15_50ns( int photonIdx );
bool isLoosePhoton_Spring15_25ns( int photonIdx );
bool isMediumPhoton_Spring15_25ns( int photonIdx );
bool isTightPhoton_Spring15_25ns( int photonIdx );
bool isLoosePhotonPOG_Spring16( int photonIdx );
bool isMediumPhotonPOG_Spring16( int photonIdx );
bool isTightPhotonPOG_Spring16( int photonIdx );
bool isLoosePhotonPOG_Fall17V2( int photonIdx );
bool isMediumPhotonPOG_Fall17V2( int photonIdx );
bool isTightPhotonPOG_Fall17V2( int photonIdx );
bool passTriggerEmu( int photonIdx );

#endif
