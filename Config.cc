#include "Config.h"
#include <iostream>

GlobalConfig gconf;

// Centrally maintained configurations
// Pass in the datasetname or filename (that contains the dataset info)
void GlobalConfig::GetConfigsFromDatasetName(TString dsname)
{
// bool isData = dsname.Contains("Run201") || dsname.Contains("run2_data");

  if (dsname.Contains("Run2016")
      || dsname.Contains("Moriond17")
      || dsname.Contains("RunIISummer16")
      || dsname.Contains("RunIISpring16")
      || dsname.Contains("run2_data2016")
      || dsname.Contains("run2_moriond17")
      ) year = 2016;
  if (dsname.Contains("Run2017")
      || dsname.Contains("RunIIFall17")
      || dsname.Contains("_mc2017_")
      || dsname.Contains("run2_mc2017")
      ) year = 2017;
  if (dsname.Contains("Run2018")
      || dsname.Contains("RunIISpring18")
      || dsname.Contains("RunIISummer18")
      || dsname.Contains("run2_mc2018")
      ) year = 2018;

  std::cout << ">>> ------------ GlobalConfig ------------" << std::endl;
  if (year <= 0) {
    std::cout << ">>> [!] Couldn't figure out year, so setting it to 2017. Make sure this is what you want!" << std::endl;
    year = 2017;
  } else {
    std::cout << ">>> Figured out that the year is " << year << "." << std::endl;
  }
  std::cout << ">>> --------------------------------------" << std::endl;

  if (dsname.Contains("80X")
      || dsname.Contains("Moriond17")
      || dsname.Contains("RunIISummer16MiniAODv2")
      || dsname.Contains("03Feb2017")  // 2016 data for Moriond17
      ) cmssw_ver = 80;
  if (dsname.Contains("94X")
      || dsname.Contains("RunIIFall17MiniAODv2")
      || dsname.Contains("31Mar2018")  // 2017 data ReReco
      || dsname.Contains("17Jul2018")  // 2016 data ReReco
      ) cmssw_ver = 94;
  if (dsname.Contains("101X")
      ) cmssw_ver = 101;

  if (year == 2016 && cmssw_ver == 80) {
    ea_version = 1;
    btag_disc_wp = 0.6324;
    jecEraB = jecEraC = jecEraD = "Summer16_23Sep2016BCDV4_DATA";
    jecEraE = jecEraF = "Summer16_23Sep2016EFV4_DATA";
    jecEraG = "Summer16_23Sep2016GV4_DATA";
    jecEraH = "Summer16_23Sep2016HV4_DATA";
    jecEraMC = "Summer16_23Sep2016V4_MC";
    jecEraFS = "Spring16_FastSimV1";

    // B-tag working points
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    // Data: 23Sep2016 ReReco for B-G datasets
    // Monte Carlo: RunIISummer16 for Fullsim and Spring16 for Fastsim

    WP_DEEPCSV_TIGHT  = 0.8958;
    WP_DEEPCSV_MEDIUM = 0.6324;
    WP_DEEPCSV_LOOSE  = 0.2219;
    fn_btagSF_DeepCSV = "DeepCSV_Moriond17_B_H.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv";

    WP_CSVv2_TIGHT  = 0.9535;
    WP_CSVv2_MEDIUM = 0.8484;
    WP_CSVv2_LOOSE  = 0.5426;
    fn_btagSF_CSVv2 = "CSVv2_Moriond17_B_H.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv";
  }
  if (year == 2016 && cmssw_ver == 94) {
    ea_version = 1;
    btag_disc_wp = 0.6324;
    jecEraB = jecEraC = jecEraD = "Summer16_07Aug2017BCD_V11_DATA";
    jecEraE = jecEraF = "Summer16_07Aug2017EF_V11_DATA";
    jecEraG = jecEraH = "Summer16_07Aug2017GH_V11_DATA";
    jecEraMC = "Summer16_07Aug2017_V11_MC";
    jecEraFS = "Spring16_FastSimV1";      // to be updated

    // B-tag working points
    // 94X WPs shall be very close to those in 80X, if not the same
    WP_DEEPCSV_TIGHT  = 0.8958;
    WP_DEEPCSV_MEDIUM = 0.6324;
    WP_DEEPCSV_LOOSE  = 0.2219;
    fn_btagSF_DeepCSV = "DeepCSV_Moriond17_B_H.csv";              // to be updated
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv"; // to be updated

    WP_CSVv2_TIGHT  = 0.9535;
    WP_CSVv2_MEDIUM = 0.8484;
    WP_CSVv2_LOOSE  = 0.5426;
    fn_btagSF_CSVv2 = "CSVv2_Moriond17_B_H.csv";               // to be updated
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv";  // to be updated
  }
  else if (year == 2017) {
    ea_version = 4;
    btag_disc_wp = 0.4941;
    jecEraB = "Fall17_17Nov2017B_V32_DATA";
    jecEraC = "Fall17_17Nov2017C_V32_DATA";
    jecEraD = jecEraE = "Fall17_17Nov2017DE_V32_DATA";
    jecEraF = "Fall17_17Nov2017F_V32_DATA";
    jecEraMC = "Fall17_17Nov2017_V32_MC";

    // B-tag working points
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    // Data: 17Nov2017 ReReco for B-F dataset
    // Monte Carlo: RunIIFall17

    WP_DEEPCSV_TIGHT  = 0.8001;
    WP_DEEPCSV_MEDIUM = 0.4941;
    WP_DEEPCSV_LOOSE  = 0.1522;
    fn_btagSF_DeepCSV = "DeepCSV_94XSF_V3_B_F.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv"; // to be updated

    WP_CSVv2_TIGHT  = 0.9693;
    WP_CSVv2_MEDIUM = 0.8838;
    WP_CSVv2_LOOSE  = 0.5803;
    fn_btagSF_CSVv2 = "CSVv2_94XSF_V2_B_F.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv"; // to be updated
  }
  else if (year == 2018) {
    // everything to be updated
    ea_version = 4;
    btag_disc_wp = 0.4184;
    jecEraA = jecEraB = jecEraC = jecEraD = "Fall17_17Nov2017F_V32_DATA"; // to be updated
    jecEraMC = "Fall17_17Nov2017_V32_MC"; // to be updated

    // B-tag working points
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    WP_DEEPCSV_TIGHT  = 0.7527;
    WP_DEEPCSV_MEDIUM = 0.4184;
    WP_DEEPCSV_LOOSE  = 0.1241;
    fn_btagSF_DeepCSV = "DeepCSV_94XSF_V3_B_F.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv"; // to be updated

    WP_CSVv2_TIGHT  = 0.9693; // FIXME I believe CSVv2 is no longer supported for 2018
    WP_CSVv2_MEDIUM = 0.8838;
    WP_CSVv2_LOOSE  = 0.5803;
    fn_btagSF_CSVv2 = "CSVv2_94XSF_V2_B_F.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv"; // to be updated
  }

}
