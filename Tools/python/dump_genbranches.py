import os
import ROOT as r
from looper_with_core import Looper

class GenWeightChecker(Looper):

    def before_loop(self):
        def select(x):
            return any(x.startswith(y) for y in ["gen_LHEweight", "genHEPMC", "genps_weight",])
        self.toplot = filter(select,looper.get_branch_names())
        self.nominal = "genps_weight"
        self.hists = {}
        for name in self.toplot:
            self.hists[name] = r.TH1F(name,"{}".format(name,self.nominal),100,0.0,3.0)

    def after_loop(self):
        os.system("mkdir -p plots")
        c1 = r.TCanvas()
        for name in self.hists:
            h = self.hists[name]
            h.Draw("histe")
            fname = "plots/h_{}.png".format(name)
            print ">>> Made {}".format(fname)
            c1.SaveAs(fname)

    def process_event(self, ievent):
        for name in self.toplot:
            val = getattr(r.cms3,name)()
            self.hists[name].Fill(val)

if __name__ == "__main__":

    ch = r.TChain("Events")
    ch.Add("/hadoop/cms/store/group/snt/run2_mc2018/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-02/merged_ntuple_1.root")
    # ch.Add("/hadoop/cms/store/group/snt/run2_mc2018//QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2_MINIAODSIM_CMS4_V10-02-02//merged_ntuple_11.root")

    looper = GenWeightChecker(ch)
    looper.loop()

