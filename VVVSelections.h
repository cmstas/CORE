#include "CMS3.h"
#include "Base.h"

#ifndef VVVSELECTIONS_H
#define VVVSELECTIONS_H


bool passElectronSelection_VVV(int index, id_level_t id_string);
bool passElectronSelection_VVV_v1(int index, id_level_t id_string, bool vetoTransition, bool eta24);

bool passMuonSelection_VVV(int index, id_level_t id_string);
bool passMuonSelection_VVV_v1(int index, id_level_t id_string, bool vetoTransition, bool eta24);


#endif
