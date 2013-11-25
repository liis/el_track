import ROOT
from tree_variables import var_list, var_type

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ptmode', dest='ptmode', choices=["Pt1","Pt10","Pt100","PtFlat"], required=True, help= "pt cut in analyzed dataset")

args = parser.parse_args()

indir = "$WORKING_DIR/tree_to_histo/input_trees/"
infile = "trackValTree_el" + args.ptmode + ".root"

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 10000
nEvt = t.GetEntries()

print "Analyzing dataset for " + args.ptmode + " with " + str(nEvt) + " events"

#### initialize tree #######
vt = dict([ (v, var_type(v)) for v in var_list ]) #associate proper data type for every variable in the tree

t.SetBranchStatus("*",0)
for v in var_list: # relate vList to the tree
    t.SetBranchStatus(v, 1)
    t.SetBranchAddress(v, vt[v])
    t.AddBranchToCache(v, ROOT.kTRUE)
    
t.StopCacheLearningPhase()

neta = 100
mineta = -2.5
maxeta = 2.5

h_reco_eta = ROOT.TH1F("nreco_eta","nrReco vs eta", neta, mineta, maxeta)
h_fake_eta = ROOT.TH1F("nfake_eta","nrFake vs eta", neta, mineta, maxeta)
h_sim_eta = ROOT.TH1F("nsim_eta", "nrSim vs eta", neta, mineta, maxeta)
h_sim_to_reco_match_eta = ROOT.TH1F("nsimtoreco_eta", "nr sim to reco matched vs eta", neta,mineta,maxeta)

h_fakerate_eta = ROOT.TH1F("fake_rate_eta","Fake rate", neta, mineta, maxeta)
h_eff_eta = ROOT.TH1F("eff_eta","Efficiency", neta, mineta, maxeta)

npt = 100
minpt = 0
maxpt = 200

h_reco_pt = ROOT.TH1F("nreco_pt","nrReco vs pt", npt, minpt, maxpt)
h_fake_pt = ROOT.TH1F("nfake_pt","nrFake vs pt", npt, minpt, maxpt)
h_sim_pt = ROOT.TH1F("nsim_pt", "nrSim vs pt", npt, minpt, maxpt)
h_sim_to_reco_match_pt = ROOT.TH1F("nsimtoreco_pt", "nr sim to reco matched vs pt", npt,minpt,maxpt)

h_fakerate_pt = ROOT.TH1F("fake_rate_pt","Fake rate", npt, minpt, maxpt)
h_eff_pt = ROOT.TH1F("eff_pt","Efficiency", npt, minpt, maxpt)

# Event loop
for i in range(nEvt):
    if i % report_every == 0: print "Event nr: " + str(i)
    t.LoadTree(i)
    t.GetEntry(i)
    for it_p in range( vt['np_reco'][0]): # loop over reco-tracks
        h_reco_eta.Fill(vt['reco_eta'][it_p])
        h_reco_pt.Fill(vt['reco_pt'][it_p])

    for it_p in range( vt['np_fake'][0]): # loop over fake traks
        h_fake_eta.Fill(vt['fake_eta'][it_p])
        h_fake_pt.Fill(vt['fake_pt'][it_p])

    for it_p in range( vt['np_gen'][0]):
        if(vt['gen_pt'][it_p] > 99 and vt['gen_pt'][it_p] < 101):
            h_sim_eta.Fill(vt['gen_eta'][it_p])
            h_sim_pt.Fill(vt['gen_pt'][it_p])

    for it_p in range( vt['np_gen_toReco'][0]):
        if(vt['gen_pt'][it_p] >99 and vt['gen_pt'][it_p] < 101):
            h_sim_to_reco_match_eta.Fill(vt['gen_matched_eta'][it_p])
            h_sim_to_reco_match_pt.Fill(vt['gen_matched_pt'][it_p])

# efficiency and fake rate wrt eta
for i in range(1,neta+1): 
    if not h_reco_eta.GetBinContent(i) == 0:
        h_fakerate_eta.SetBinContent(i, h_fake_eta.GetBinContent(i)/h_reco_eta.GetBinContent(i) )
    else:
        h_fakerate_eta.SetBinContent(i,0)
        
    if not h_sim_eta.GetBinContent(i) == 0:
        h_eff_eta.SetBinContent(i, h_sim_to_reco_match_eta.GetBinContent(i)/h_sim_eta.GetBinContent(i) )
    else:
        h_sim_eta.SetBinContent(i,0)

# efficiency and fake rate wrt pt
for i in range(1,npt+1):
    if not h_reco_pt.GetBinContent(i) == 0:
        h_fakerate_pt.SetBinContent(i, h_fake_pt.GetBinContent(i)/h_reco_pt.GetBinContent(i) )
    else:
        h_fakerate_pt.SetBinContent(i,0)
        
    if not h_sim_pt.GetBinContent(i) == 0:
        h_eff_pt.SetBinContent(i, h_sim_to_reco_match_pt.GetBinContent(i)/h_sim_pt.GetBinContent(i) )
    else:
        h_sim_pt.SetBinContent(i,0)
                                                    

outfile = "histograms/trackValHistograms" + args.ptmode  + ".root"    
o = ROOT.TFile(outfile,"recreate")
h_reco_eta.Write()
h_fake_eta.Write()
h_sim_eta.Write()
h_sim_to_reco_match_eta.Write()

h_reco_pt.Write()
h_fake_pt.Write()
h_sim_pt.Write()
h_sim_to_reco_match_pt.Write()

h_fakerate_eta.Write()
h_eff_eta.Write()
h_fakerate_pt.Write()
h_eff_pt.Write()

o.Close()
