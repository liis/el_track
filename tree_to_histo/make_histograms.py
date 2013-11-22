import ROOT
from tree_variables import var_list, var_type

indir = "$WORKING_DIR/tree_to_histo/input_trees/"
infile = "trackValTree_elPt100.root"

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 1000

nEvt = t.GetEntries()

print "Analyzing " + str(nEvt) + " events"

#### initialize tree #######
t.SetBranchStatus("*",0)

vList = dict([ (v, var_type(v)) for v in var_list ]) #associate proper data type for every variable in the tree

for v in var_list: # relate vList to the tree
    t.SetBranchStatus(v, 1)
    t.SetBranchAddress(v, vList[v])
    t.AddBranchToCache(v, ROOT.kTRUE)
    
t.StopCacheLearningPhase()

# Event loop
for i in range(nEvt):
    if i % report_every == 0: print "Event nr: " + str(i)
    t.LoadTree(i)
    t.GetEntry(i)
    print vList['np_reco'][0]
        
