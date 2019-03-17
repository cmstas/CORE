import time
import os
import sys
import numpy as np

# Want to be able to pass numpy arrays of jet information and get back
# a numpy array of JECs, b-tag SFs, ...

np.set_printoptions(linewidth=220,precision=3)

if __name__ == "__main__":

    # Import ROOT and load CORE .so file after it has been built with the JetCorrector stuff
    # (e.g., CMS3_CORE.so, or maybe your babymaker)
    import ROOT
    ROOT.gROOT.ProcessLine(".L ../../../liblooperBatch.so")
    # Load vectorized JEC functions
    ROOT.gROOT.ProcessLine(".L PyCORE.h")

    # Instantiate corrector with L1,L2,L3 files
    from ROOT import PyCOREBridge
    pcb = PyCOREBridge()
    pcb.AddL1JECFile("../../Tools/jetcorr/data/run2_25ns/Fall17_FastsimV1/Fall17_FastsimV1_L1FastJet_AK4PFchs.txt")
    pcb.AddL2JECFile("../../Tools/jetcorr/data/run2_25ns/Fall17_FastsimV1/Fall17_FastsimV1_L2Relative_AK4PFchs.txt")
    pcb.AddL3JECFile("../../Tools/jetcorr/data/run2_25ns/Fall17_FastsimV1/Fall17_FastsimV1_L3Absolute_AK4PFchs.txt")
    pcb.MakeJetCorrector()

    # Make `N` randomly distributed jets
    N = 100000
    pts = 50.+100*np.random.random(N)
    etas = -2.4+4.8*np.random.random(N)
    rhos = np.clip(np.random.normal(25,5,N),10,40)
    areas = 0.4*np.ones(N)


    # Make an array to fill with correction values, and then get the JECs
    t0 = time.time()
    outs = np.ones(N)
    pcb.GetJECs(len(pts),pts,etas,rhos,areas,outs)
    t1 = time.time()
    print outs

    # Should be ~0.4MHz
    print "Calculated {} JEC values in {:.2f}s ({:.2f}MHz)".format(len(outs),t1-t0,1.e-6*len(outs)/(t1-t0))
