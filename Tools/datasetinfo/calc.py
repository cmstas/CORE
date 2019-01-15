import os
import json
import glob
import sys
import datetime

def get_datasettag(metadata):
    with open(metadata,"r") as fhin:
        data = json.load(fhin)
        tag = data["tag"]
        dataset = data["dataset"]
        return (dataset,tag)

def get_line(metadata):
    if not metadata.endswith("metadata.json"): metadata += "/metadata.json"
    with open(metadata,"r") as fhin:
        data = json.load(fhin)
        ijob_to_nevents = data["ijob_to_nevents"]
        tag = data["tag"]
        dataset = data["dataset"]
        kfact = data["kfact"]
        xsec = data["xsec"]
        efact = data["efact"]
        thedir = metadata.rsplit("/",1)[0]
        file_indices = map(lambda x: x.rsplit("_",1)[-1].split(".")[0], glob.glob(thedir+"/*.root"))
        nevts_total = 0
        nevts_eff_total = 0
        good_pairs = []
        for fi in file_indices:
            try:
                nevts, nevts_eff = ijob_to_nevents[fi]
            except:
                print "ERROR with fi={}".format(fi)
            good_pairs.append([fi,(nevts,nevts_eff)])
            nevts_total += nevts
            nevts_eff_total += nevts_eff
        xsec_total = kfact*xsec*efact
        scale1fb = 1000.0*xsec_total/nevts_eff_total

        info = ",".join(["{}|{}|{}".format(idx,nevt,(nevt-nevteff)/2) for idx,(nevt,nevteff) in good_pairs])

        return "{:155s}  {:17s}  {:9}  {:9}  {:8.5g}  {:10.5g} {}".format(dataset,tag,nevts_total,nevts_eff_total,xsec_total,scale1fb,info)

if __name__ == "__main__":

    sampledirs = []

    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc2017/*09-04-17*/")
    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc2017/*09-04-19*/")
    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc2017/*09-04-20*/")
    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc2018/*10-02-02*/")
    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc2016_94x/*09-04-17*/")
    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc2016_94x/*09-04-17*/")
    # sampledirs += glob.glob("/hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/*CMS4_V00-00-02_2017Sep27*/")
    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc2016_cms4/*/")
    sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc*/*10-02-05/")
    sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc*/*10-02-04/")
    # sampledirs += glob.glob("/hadoop/cms/store/group/snt/run2_mc*/*10-02-04/")

    # print get_line("/hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2_MINIAODSIM_CMS4_V00-00-02_2017Sep27/metadata.json")
    # print get_line("/hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2_MINIAODSIM_CMS4_V00-00-02_2017Sep27/metadata.json")
    # print get_line("/hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_MINIAODSIM_CMS4_V00-00-02_2017Sep27/metadata.json")
    # print get_line("/hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_MINIAODSIM_CMS4_V00-00-02_2017Sep27/metadata.json")
    # print get_line("/hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/tZq_ll_4f_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_MINIAODSIM_CMS4_V00-00-02_2017Sep27/metadata.json")
    # print get_line("/hadoop/cms/store/group/snt/run2_mc2017//TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V09-04-19/")
    # print get_line("/hadoop/cms/store/group/snt/run2_mc2017/TTTW_TuneCP5_13TeV-madgraph-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V09-04-19//metadata.json")
    # sampledirs = glob.glob("/hadoop/cms/store/group/snt/run2_mc2018/*V10-02-00*/")
    # blah

    alreadydone = set()
    with open("scale1fbs.txt", "r") as fhin:
        for line in fhin:
            parts = line.strip().split()
            if len(parts) < 6: continue
            dsname = parts[0]
            tag = parts[1]
            alreadydone.add((dsname,tag))


    print "# file created on %s" % (datetime.datetime.now())
    print "# {:153s}  {:17s}  {:9s}  {:9s}  {:8s}  {:10s}".format("dataset","tag","nevts_tot","nevts_eff","xsec","scale1fb")
    for sampledir in sampledirs:
        metadata = sampledir + "/metadata.json"
        if not os.path.exists(metadata): continue
        try:
            if get_datasettag(metadata) in alreadydone: continue
            print get_line(metadata)
        except:
            sys.stderr.write("[!] Error with {}\n".format(metadata))
