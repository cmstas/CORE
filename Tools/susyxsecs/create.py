import numpy as np
import ROOT as r

# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVgluglu
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVsquarkantisquark
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVstopsbottom
# Copy paste tables from the twikis into text files. The format should look like
# 500. GeV	0.338E+02 ± 7.86 %
# 505. GeV	0.319E+02 ± 7.89 %
# 510. GeV	0.301E+02 ± 7.92 %

def parse(x):
    with open(x,"r") as fh:
        data = []
        for line in fh.readlines():
            if not line.strip(): continue
            parts = line.split()
            mass = float(parts[0])
            xsec = float(parts[2])
            pcterr = float(parts[4])
            xsecerr = pcterr/100. * xsec
            data.append([mass,xsec,xsecerr])
    return np.array(data)

def make_hist(data, name):
    binwidth = data[:,0][1]-data[:,0][0]
    print name,name,len(data),data[:,0].min()-0.5*binwidth,data[:,0].max()+0.5*binwidth
    h = r.TH1F(name,name,len(data),data[:,0].min()-0.5*binwidth,data[:,0].max()+0.5*binwidth)
    for mass,xsec,xsecerr in data:
        ib = h.FindBin(mass)
        h.SetBinContent(ib,xsec)
        h.SetBinError(ib,xsecerr)
    return h

f = r.TFile("xsec_susy_13tev_run2.root","RECREATE")
h_xsec_gluino = make_hist(parse("gluino.txt"), "h_xsec_gluino")
h_xsec_gluino.Write()
h_xsec_squark = make_hist(parse("squark.txt"), "h_xsec_squark")
h_xsec_squark.Write()
h_xsec_stop = make_hist(parse("stop.txt"), "h_xsec_stop")
h_xsec_stop.Write()
f.Close()

