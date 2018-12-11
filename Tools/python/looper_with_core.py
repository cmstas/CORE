import os
import ROOT as r
import time
from tqdm import tqdm

class Looper(object):
    """
    Import, subclass, override before_loop, after_loop, and process_event, then
    >>> l = Looper(mytchain)
    >>> l.loop()
    """
    def __init__(self, chain, extra_includes=[]):
        self.chain = chain
        self.coredir = "../../"
        if not os.path.exists("{}/CMS3_CORE.so".format(self.coredir)):
            raise Exception("Hey, you need to compile in order to get CMS3_CORE.so")
        for include in ["CMS3.cc", "JetSelections.h","MetSelections.h","ElectronSelections.h"]+extra_includes:
            r.gInterpreter.ProcessLine('#include "{}/{}"'.format(self.coredir,include))
        r.gSystem.Load("{}/CMS3_CORE.so".format(self.coredir))

    def get_branch_names(self):
        symbols = set(filter(lambda x:not x.startswith("_"), dir(r.cms3)))
        notbranches = set(["GetEntry","Init","LoadAllBranches","progress"])
        return sorted(list(symbols-notbranches))

    def before_loop(self):
        """
        Subclass and override this
        """
        self.h1 = r.TH1F("h1","my hist", 10, 0., 500.)
        self.h2 = r.TH1F("h2","my hist", 10, 0., 1.0)

    def after_loop(self):
        """
        Subclass and override this
        """
        print list(self.h1)
        print list(self.h2)

    def process_event(self, ievent):
        """
        Subclass and override this
        """
        # Get a branch directly
        self.h1.Fill(r.cms3.evt_pfmet())
        # Use a CORE function
        triplet = r.getPrefireInfo(2017) # std::tuple
        self.h2.Fill(triplet._0)

    def loop(self, progress=True):
        r.cms3.Init(self.chain)
        self.before_loop()
        if not progress:
            iterable = range(self.chain.GetEntries())
        else:
            iterable = tqdm(range(self.chain.GetEntries()),total=self.chain.GetEntries())
        for ievent in iterable:
            r.cms3.GetEntry(ievent)
            self.process_event(ievent)
        self.after_loop()


if __name__ == "__main__":

    ch = r.TChain("Events")
    ch.Add("/hadoop/cms/store/group/snt/run2_mc2018/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-02/merged_ntuple_1.root")

    looper = Looper(ch)
    looper.loop()


