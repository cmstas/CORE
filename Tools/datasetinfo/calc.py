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
        for fi in file_indices:
            try:
                nevts, nevts_eff = ijob_to_nevents[fi]
            except:
                print "ERROR with fi={}".format(fi)
            nevts_total += nevts
            nevts_eff_total += nevts_eff
        xsec_total = kfact*xsec*efact
        scale1fb = 1000.0*xsec_total/nevts_eff_total

        info = ",".join(["{}|{}|{}".format(idx,nevt,(nevt-nevteff)/2) for idx,(nevt,nevteff) in data["ijob_to_nevents"].items()])

        return "{:155s}  {:17s}  {:9}  {:9}  {:8.5g}  {:10.5g} {}".format(dataset,tag,nevts_total,nevts_eff_total,xsec_total,scale1fb,info)

if __name__ == "__main__":

    sampledirs = glob.glob("/hadoop/cms/store/group/snt/run2_mc2017/*09-04-13*/")

    alreadydone = set()
    with open("scale1fbs.txt", "r") as fhin:
        for line in fhin:
            parts = line.strip().split()
            if len(parts) != 6: continue
            dsname = parts[0]
            tag = parts[1]
            alreadydone.add((dsname,tag))


    print "# file created on %s" % (datetime.datetime.now())
    print "# {:153s}  {:17s}  {:9s}  {:9s}  {:8s}  {:10s}".format("dataset","tag","nevts_tot","nevts_eff","xsec","scale1fb")
    for sampledir in sampledirs:
        metadata = sampledir + "/metadata.json"
        if not os.path.exists(metadata): continue
        if get_datasettag(metadata) in alreadydone: continue
        print get_line(metadata)
