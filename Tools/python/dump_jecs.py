import os
import sys
import numpy as np

np.set_printoptions(linewidth=220,precision=3)

pts = np.array([10.,15.,20.,25.,30.,35.,40.,45.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200,240,280,320,360,400,500])
etas = np.array([0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5])
etas = np.concatenate([-etas[::-1],etas])

makeJetCorrector, r = None, None
def load():
    global makeJetCorrector, r
    import ROOT as r
    r.gROOT.ProcessLine(".L ../../../CORE/CMS3_CORE.so")
    for include in [
            "CMS3.cc",
            "Tools/JetCorrector.h",
            "Tools/jetcorr/JetCorrectionUncertainty.h",
            "Tools/jetcorr/SimpleJetCorrectionUncertainty.h",
            ]:
        r.gInterpreter.ProcessLine('#include "../../../CORE/%s"' % include)
    from ROOT import makeJetCorrector

def dump():
    if not makeJetCorrector:
        load()

    fnames = r.std.vector("string")()
    fnames.push_back("../../../CORE/Tools/jetcorr/data/run2_25ns/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK4PFchs.txt")
    fnames.push_back("../../../CORE/Tools/jetcorr/data/run2_25ns/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK4PFchs.txt")
    fnames.push_back("../../../CORE/Tools/jetcorr/data/run2_25ns/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK4PFchs.txt")

    jc = makeJetCorrector(fnames)

    data = []
    for pt in pts:
        for eta in etas:
            jc.setRho(22.0)
            jc.setJetA(0.5)
            jc.setJetPt(pt)
            jc.setJetEta(eta)
            corrs = jc.getSubCorrections()
            # last value is final correction
            # pt, eta, L1, L2L3, L1L2L3
            data.append([pt,eta,corrs[0],corrs[-1]/corrs[0],corrs[-1]])
    data = np.array(data)
    return data


if __name__ == "__main__":

    data = dump()
    print np.histogram2d(data[:,0],data[:,1], weights=data[:,2],bins=[pts,etas])
