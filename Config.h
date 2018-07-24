#ifndef CONFIG_H
#define CONFIG_H
#include <string>

#include "Tools/JetCorrector.h"

/*
   Need to check the year in a CORE function (or your own)? `gconf` is a global instance.
       if (gconf.year == 2017) { ... }
   Want to update the year in your babymaker?
       gconf.year = 2018;
   Need another variable? Add a line with a dummy default. If an analysis uses
   a variable from this config, it's up to them to make sure the variables are set properly
   at runtime. Leaving dummy values here forces you to do that.
*/

class GlobalConfig {
    public:
        unsigned int year = 0;
        unsigned int cmssw_ver = 0; // 74, 80, 94, 101, ...
        std::string analysis = ""; 
        float btag_disc_wp = -1;
        int ea_version = -1;

        // MultiIso WPs for SS
        float multiiso_el_minireliso = -1;
        float multiiso_el_ptratio = -1;
        float multiiso_el_ptrel = -1;
        float multiiso_mu_minireliso = -1;
        float multiiso_mu_ptratio = -1;
        float multiiso_mu_ptrel = -1;

        // JECs
        FactorizedJetCorrector * jet_corrector_L1 = 0;
        FactorizedJetCorrector * jet_corrector_L2L3 = 0;
        FactorizedJetCorrector * jet_corrector_L1L2L3 = 0;

        // ...
};

#ifndef __CINT__
extern GlobalConfig gconf;
#endif

#endif
