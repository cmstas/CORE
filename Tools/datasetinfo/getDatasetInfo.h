#ifndef CORE_TOOLS_GETDATASETINFO_H
#define CORE_TOOLS_GETDATASETINFO_H

#include <string>
#include <unordered_map>

using namespace std;
using namentag = std::pair<string,string>;

class DatasetInfoFromFile {
 public:
  struct datasetInfo {
    unsigned int nevts_tot;
    unsigned int nevts_eff;
    float xsec;
    float scale1fb;
  };
  DatasetInfoFromFile();  // This will load the default datset list defined in defaultList.icc
  DatasetInfoFromFile(const string filename) { loadFromFile(filename); }
  void loadFromFile(const string filename);
  void checkEntryExist(const string datasetname, const string cmstag);
  float getScale1fbFromFile(const string datasetname, const string cmstag);
  float getXsecFromFile(const string datasetname, const string cmstag);
  unsigned int getnEventsTotalFromFile(const string datasetname, const string cmstag);
  unsigned int getnEventsEffectiveFromFile(const string datasetname, const string cmstag);
  bool isEmpty() { return dslist_.empty(); };
  size_t numberOfEntries() { return dslist_.size(); };
  void update(DatasetInfoFromFile::datasetInfo* df, string info, string dsname, string tag, bool verbose=false);

 private:
  std::unordered_map<string, datasetInfo> dslist_;
};

#endif  // GETDATASETINFO_H
