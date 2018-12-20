#ifndef CORE_TOOLS_GETDATASETINFO_H
#define CORE_TOOLS_GETDATASETINFO_H

#include <string>
#include <unordered_map>


class DatasetInfoFromFile {
 public:
  struct datasetInfo {
    unsigned int nevts_tot;
    unsigned int nevts_eff;
    float xsec;
    float scale1fb;
  };
  DatasetInfoFromFile();  // This will load the default dataset list defined in defaultList.icc
  DatasetInfoFromFile(const std::string filename) { loadFromFile(filename); }
  void loadFromFile(const std::string filename, bool verbose=false);
  void checkEntryExist(const std::string datasetname, const std::string cmstag);
  bool doesEntryExist(const std::string datasetname, const std::string cmstag);
  float getScale1fbFromFile(const std::string datasetname, const std::string cmstag);
  float getXsecFromFile(const std::string datasetname, const std::string cmstag);
  unsigned int getnEventsTotalFromFile(const std::string datasetname, const std::string cmstag);
  unsigned int getnEventsEffectiveFromFile(const std::string datasetname, const std::string cmstag);
  bool isEmpty() { return dslist_.empty(); };
  size_t numberOfEntries() { return dslist_.size(); };
  void update(DatasetInfoFromFile::datasetInfo* df, std::string info, std::string dsname, std::string tag, bool verbose);

 private:
  std::unordered_map<std::string, datasetInfo> dslist_;
};

#endif  // GETDATASETINFO_H
