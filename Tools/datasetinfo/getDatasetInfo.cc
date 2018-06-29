#include <iostream>
#include <fstream>
#include <sstream>
#include "getDatasetInfo.h"
#include "defaultList.icc"

using namespace std;

DatasetInfoFromFile::DatasetInfoFromFile() : dslist_(default_dslist) {}

void DatasetInfoFromFile::update(DatasetInfoFromFile::datasetInfo* df, string info, string dsname, string tag, bool verbose) {
    istringstream iss(info);
    string token;
    int nevents_total = 0;
    int nevents_negative = 0;
    while (getline(iss, token, ',')) {
        istringstream parts(token);
        string part;
        int ipart = 0;
        int idx = -1;
        while (getline(parts, part, '|')) {
            switch (ipart) {
                case 0: idx = stoi(part); break;
                case 1: nevents_total += stoi(part); break;
                case 2: nevents_negative += stoi(part); break;
            }
            ipart += 1;
        }
    }
    int old_nevts_tot = df->nevts_tot;
    int new_nevts_tot = nevents_total;
    int old_nevts_eff = df->nevts_eff;
    int new_nevts_eff = nevents_total-2*nevents_negative; // Ntot = N(+) + N(-), Neff = N(+) - N(-)
    float old_scale1fb = df->scale1fb;
    float new_scale1fb = 1000. * df->xsec / new_nevts_eff;
    if (verbose && ((old_nevts_tot != new_nevts_tot) || (old_nevts_eff != new_nevts_eff))) {
        std::cout << "=== " << dsname << " [" << tag << "] ===" << std::endl;
        printf("nevents total: %i -> %i [%+i]\n", old_nevts_tot, new_nevts_tot, new_nevts_tot-old_nevts_tot);
        printf("nevents effective: %i -> %i [%+i]\n", old_nevts_eff, new_nevts_eff, new_nevts_eff-old_nevts_eff);
        printf("scale1fb: %g -> %g [%+.2f]\n", old_scale1fb, new_scale1fb, 100.*(new_scale1fb-old_scale1fb)/old_scale1fb);
        std::cout << "========================" << std::endl;
    }
    df->nevts_tot = new_nevts_tot;
    df->nevts_eff = new_nevts_eff;
    df->scale1fb = new_scale1fb;
}

void DatasetInfoFromFile::loadFromFile(const string filename) {
  ifstream ifile(filename);
  if (!ifile) throw std::invalid_argument("Dataset info file " + filename + " does not exist!");
  string dsname, cmstag;
  datasetInfo dsinfo;
  string line;
  while (getline(ifile, line)) {
    if (line.find("# dataset") == 0) {}    // check the ordering of the variables <-- implement later
    if (line[0] == '#') continue;
    istringstream lstream(line);
    string entry;
    for (int i = 0; lstream >> std::ws >> entry; ++i) {
      switch (i) {
      case 0: dsname = entry; break;
      case 1: cmstag = entry; break;
      case 2: dsinfo.nevts_tot = stoul(entry); break;
      case 3: dsinfo.nevts_eff = stoul(entry); break;
      case 4: dsinfo.xsec = stof(entry); break;
      case 5: dsinfo.scale1fb = stof(entry); break;
      case 6: update(&dsinfo, entry, dsname, cmstag, true); break;
      case 7: throw std::invalid_argument("Dataset info file has too many arguments! Please update this code!");
      }
    }
    dslist_[cmstag+dsname] = dsinfo; // insert info into the dataset list
  }
}

void DatasetInfoFromFile::checkEntryExist(const string dsname, const string cmstag) {
  if (isEmpty())
    throw std::invalid_argument("The dataset list is empty. Please check if the dataset info file is properly loaded.");
  if (dslist_.find(cmstag+dsname) == dslist_.end())
    throw std::invalid_argument("The dataset " + dsname + " with cmstag: " + cmstag + " cannot be found in the list. Please update the dataset info file.");
}

float DatasetInfoFromFile::getScale1fbFromFile(const string dsname, const string cmstag) {
  checkEntryExist(dsname, cmstag);
  return dslist_[cmstag+dsname].scale1fb;
}

float DatasetInfoFromFile::getXsecFromFile(const string dsname, const string cmstag) {
  checkEntryExist(dsname, cmstag);
  return dslist_[cmstag+dsname].xsec;
}

unsigned DatasetInfoFromFile::getnEventsTotalFromFile(const string dsname, const string cmstag) {
  checkEntryExist(dsname, cmstag);
  return dslist_[cmstag+dsname].nevts_tot;
}

unsigned DatasetInfoFromFile::getnEventsEffectiveFromFile(const string dsname, const string cmstag) {
  checkEntryExist(dsname, cmstag);
  return dslist_[cmstag+dsname].nevts_eff;
}
