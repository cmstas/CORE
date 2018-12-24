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
  std::unordered_map<std::string, DatasetInfoFromFile::datasetInfo>::const_iterator checkEntryExist(const std::string datasetname, const std::string cmstag) const;
  bool doesEntryExist(const std::string datasetname, const std::string cmstag) const;
  float getScale1fbFromFile(const std::string datasetname, const std::string cmstag) const;
  float getXsecFromFile(const std::string datasetname, const std::string cmstag) const;
  unsigned int getnEventsTotalFromFile(const std::string datasetname, const std::string cmstag) const;
  unsigned int getnEventsEffectiveFromFile(const std::string datasetname, const std::string cmstag) const;
  bool isEmpty() const{ return dslist_.empty(); }
  size_t numberOfEntries() const{ return dslist_.size(); }
  void update(DatasetInfoFromFile::datasetInfo* df, std::string info, std::string dsname, std::string tag, bool verbose);
  std::unordered_map<std::string, datasetInfo> const& get_dslist(){ return dslist_; }

 private:
  std::unordered_map<std::string, datasetInfo> dslist_;
};

#endif  // GETDATASETINFO_H
