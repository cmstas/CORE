import ROOT as r
r.gROOT.ProcessLine(".L getDatasetInfo.cc")

df = r.DatasetInfoFromFile()
df.loadFromFile("scale1fbs.txt")

dsname = "/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
tag = "CMS4_V09-04-13"
ntot = df.getnEventsTotalFromFile(dsname,tag)
neff = df.getnEventsEffectiveFromFile(dsname,tag)
scale1fb = df.getScale1fbFromFile(dsname,tag)
print ntot, neff, scale1fb
