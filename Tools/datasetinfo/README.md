### scale1fb computation

scale1fb no longer lives in ntuples since there is no postprocessing step for CMS4. Instead,
metadata files containing event counts for each dataset are used to compute and store scale1fb in
`scale1fbs.txt`. A dataset name and CMStag are the only two things needed to retrieve a unique
scale1fb.

#### Updating the file
```
# update calc.py with the appropriate pattern of hadoop paths to consider
python calc.py >> scale1fbs.txt
```

#### C++ usage

They can be read in C++ via
```cpp
#include "CORE/Tools/datasetinfo/getDatasetInfo.h"

DatasetInfoFromFile df;

// Load scale1fbs/xsecs from file once
df.loadFromFile("CORE/Tools/datasetinfo/scale1fbs.txt");

// Each event, get the sign of the MC weight and multiply it by the magnitude
of the scale1fb in the text file. A scale1fb is uniquely specified by a pair of (dataset name, cms3tag), which are both embedded in the ntuples already
float sgnMCweight = ((tas::genps_weight() > 0) - (tas::genps_weight() < 0));
float scale1fb = sgnMCweight*df.getScale1fbFromFile(tas::evt_dataset()[0].Data(),tas::evt_CMS3tag()[0].Data());
float xsec = sgnMCweight*df.getXsecFromFile(tas::evt_dataset()[0].Data(),tas::evt_CMS3tag()[0].Data());
```

#### Python usage

See `test_macro.py` for example.

#### Dealing with corrupt files

Newer samples end up having a final column in the `scale1fbs.txt` file. It is a single string containing
triplets of `ifile1|nevents1|neventsnegative1,ifile2|nevents2|neventsnegative2,...`. The advantage of this is that if 
a file becomes corrupted on hadoop, only 3 things need to be done
1. Delete the file on hadoop
2. Remove the triplet for that `merged_ntuple_<i>.root` file in the `scale1fbs.txt` file
3. Commit the `scale1fbs.txt` file
The datasetinfo object will then automatically recompute the scale1fb from this reduced set of events in your looper!
Note that if this string is not modified at all, the recomputed scale1fb will match the 6th column anyway:
```
/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM                                    CMS4_V09-04-13      10881896    6820606    4.4297  0.00064946 42|189769|35328,24|286799|53367,25|236799|44265,26|283319|52730,27|260711|48563,20|254058|47140,21|234771|43677,22|283031|52894,23|278961|51901,28|243787|45864,29|235970|44431,40|259305|48366,41|283532|52662,1|235387|43943,3|236154|43881,2|231068|43144,5|282217|52818,4|283956|53088,7|283035|52869,6|282587|52969,9|284084|53015,8|258112|48248,39|276942|51664,38|259550|48581,11|260006|48566,10|236178|44114,13|275112|51383,12|236595|44153,15|282275|52393,14|260084|48647,17|282712|52606,16|210811|39424,19|283325|52745,18|260350|48675,31|261312|48593,30|261070|48748,37|283039|52932,36|236957|44504,35|253325|47098,34|258107|48085,33|259869|48523,32|236865|44048
```
