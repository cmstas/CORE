import os
import sys
import time

# pip install --user uproot
import uproot
from uproot_methods import TLorentzVectorArray

import numpy as np
# pip install --user matplotlib
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

t0 = time.time()
# Open up the file and get the Events tree
f = uproot.open("/hadoop/cms/store/group/snt/run2_mc2018//DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04/merged_ntuple_10.root")
t = f["Events"]
print("{:.2f}s to open the file and get the tree".format(time.time()-t0))

# Unfortunately, branch aliases complicate things, so make a look up table (alias->real branch name)
alias_lut = dict([(x._fName,x._fTitle) for x in t._fAliases])

def get_p4(name):
    # Given an alias like "els_p4", return a jagged array of LorentzVectors
    vx,vy,vz,vt = t.arrays([
        alias_lut[name]+".fCoordinates.fX",
        alias_lut[name]+".fCoordinates.fY",
        alias_lut[name]+".fCoordinates.fZ",
        alias_lut[name]+".fCoordinates.fT",
        ], outputtype=tuple)
    return TLorentzVectorArray.from_cartesian(vx,vy,vz,vt)

t0 = time.time()
# Get a jagged array of electrons, and include charge information
electrons = get_p4("els_p4")
electrons["charge"] = t["els_charge"].array()
print("{:.2f}s to read electron branches".format(time.time()-t0))

t0 = time.time()
# Keep only electrons with pT>20 (and note that the charges will automatically follow!)
electrons = electrons[electrons.pt>20]

# Get pairs of distinct electrons, and charges too
ee = electrons.distincts()
qee = electrons["charge"].distincts()
# Sum up the elements in the pair (i0, i1), then calculate the mass
mee = (ee.i0+ee.i1).mass
# Product of paired charges should be negative for OS pairs
isos = (qee.i0*qee.i1)<0
print("{:.2f}s to compute stuff".format(time.time()-t0))

fig,ax = plt.subplots()
# Plot invariant masses for OS electron pairs
ax.hist(mee[isos].content,bins=np.linspace(0,150,50))
fig.savefig("plot.png")
os.system("which ic && ic plot.png")

