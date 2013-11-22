import ROOT
from array import *

indir = "$WORKING_DIR/tree_to_histo/input_trees/"
infile = "trackValTree_elPt100.root"

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 1000

nEvt = t.GetEntries()

print "Analyzing " + str(nEvt) + " events"

allvars = ["np_reco", "reco_pt", "reco_eta"]

t.SetBranchStatus("*",0)
#cacheSize=10000000
#t.SetCacheSize(cacheSize)

for v in allvars:
    t.SetBranchStatus(v, 1)

vList = {"np_reco": array('i',[0]), "reco_pt": array('i',[0]*100), "reco_eta": array('i',[0]*100) }

for v in allvars:
    t.SetBranchAddress(v, vList[v])
    t.AddBranchToCache(v, ROOT.kTRUE)
    
t.StopCacheLearningPhase()

t.Draw("np_reco")

for i in range(nEvt):
    if i % report_every == 0: print "Event nr: " + str(i)
    t.LoadTree(i)
    t.GetEntry(i)
    print vList['np_reco'][0]
        
