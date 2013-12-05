import ROOT
from array import array
from tree_variables import var_list, var_type
from histlib import fill_hist_ratio, log_binning

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ptmode', dest='ptmode', choices=["Pt1","Pt10","Pt100","FlatPt"], required=True, help= "pt cut in analyzed dataset")

args = parser.parse_args()

indir = "$WORKING_DIR/tree_to_histo/input_trees/"
#infile = "trackValTree_el" + args.ptmode + ".root"
infile = "trackValTree.root" ## CLEAN UP!

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 1
nEvt = t.GetEntries()

print "Analyzing dataset for " + args.ptmode + " with " + str(nEvt) + " events"

#### initialize tree #######
vt = dict([ (v, var_type(v)) for v in var_list ]) #associate proper data type for variables in the tree

t.SetBranchStatus("*",0)
for v in var_list: # relate vList to the tree
    t.SetBranchStatus(v, 1)
    t.SetBranchAddress(v, vt[v])
    t.AddBranchToCache(v, ROOT.kTRUE)
    
t.StopCacheLearningPhase()

neta = 100
mineta = -2.5
maxeta = 2.5

h_reco_eta = ROOT.TH1F("reco_eta","nrReco vs eta", neta, mineta, maxeta)
h_fake_eta = ROOT.TH1F("fake_eta","nrFake vs eta", neta, mineta, maxeta)
h_sim_eta = ROOT.TH1F("sim_eta", "nrSim vs eta", neta, mineta, maxeta)
h_sim_to_reco_match_eta = ROOT.TH1F("nsimtoreco_eta", "nr sim to reco matched vs eta", neta,mineta,maxeta)
h_sim_passhit3_eta = ROOT.TH1F("sim_passhit3_eta", "hit eff at layer3",neta,mineta,maxeta)
h_sim_passlast_eta = ROOT.TH1F("sim_passlast_eta", "hit eff at last layer", neta, mineta, maxeta)

npt = 50
minpt = 0.1
maxpt = 200

xbins = log_binning(npt,minpt,maxpt) #log binning for pt histograms

h_reco_pt_barrel = ROOT.TH1F("nreco_pt_barrel","nrReco vs pt barrel", npt, array('d',xbins) )
h_reco_pt_endcap = ROOT.TH1F("nreco_pt_endcap","nrReco vs pt endcap", npt, array('d',xbins) )
h_reco_pt_trans = ROOT.TH1F("nreco_pt_trans","nrReco vs pt trans", npt, array('d',xbins) )

h_fake_pt_barrel = ROOT.TH1F("nfake_pt_barrel","nrFake vs pt barrel", npt, array('d', xbins) )
h_fake_pt_endcap = ROOT.TH1F("nfake_pt_endcap","nrFake vs pt endcap", npt, array('d', xbins) )
h_fake_pt_trans = ROOT.TH1F("nfake_pt_trans","nrFake vs pt trans", npt, array('d',xbins) )

h_sim_pt_barrel = ROOT.TH1F("sim_pt_barrel", "nrSim vs pt", npt, array('d',xbins))
h_sim_pt_endcap = ROOT.TH1F("sim_pt_endcap", "nrSim vs pt", npt, array('d',xbins))
h_sim_pt_trans = ROOT.TH1F("sim_pt_trans", "nrSim vs pt", npt, array('d', xbins))



h_sim_to_reco_match_pt_barrel = ROOT.TH1F("nsimtoreco_pt_barrel", "nr sim to reco matched vs pt", npt,array('d', xbins))
h_sim_to_reco_match_pt_endcap = ROOT.TH1F("nsimtoreco_pt_endcap", "nr sim to reco matched vs pt", npt, array('d',xbins))
h_sim_to_reco_match_pt_trans = ROOT.TH1F("nsimtoreco_pt_trans", "nr sim to reco matched vs pt", npt, array('d', xbins))

# Event loop
for i in range(nEvt):
    if i % report_every == 0: print "Event nr: " + str(i)
    t.LoadTree(i)
    t.GetEntry(i)
    for it_p in range( vt['np_reco'][0]): # loop over reco-tracks
        h_reco_eta.Fill(vt['reco_eta'][it_p])
        
        if abs(vt['reco_eta'][it_p]) < 0.9:
            h_reco_pt_barrel.Fill(vt['reco_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 1.4:
            h_reco_pt_trans.Fill(vt['reco_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 2.5:
            h_reco_pt_endcap.Fill(vt['reco_pt'][it_p])

    for it_p in range( vt['np_fake'][0]): # loop over fake traks
        h_fake_eta.Fill(vt['fake_eta'][it_p])

        if abs(vt['reco_eta'][it_p]) < 0.9:
            h_fake_pt_barrel.Fill(vt['fake_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 1.4:
            h_fake_pt_trans.Fill(vt['fake_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 2.5:
            h_fake_pt_endcap.Fill(vt['fake_pt'][it_p])

    for it_p in range( vt['np_gen'][0]): #loop over simulated tracks
        print "Found generated track"
        for nrhit in range(0,30):
            print vt["gen_passhit075"][it_p][nrhit]
        
        h_sim_eta.Fill(vt['gen_eta'][it_p])
        
        if(vt['gen_passhit3_eta'][it_p] > -100):
            h_sim_passhit3_eta.Fill(vt['gen_eta'][it_p])
          #  h_sim_passhit3_pt.Fill(vt['gen_pt'][it_p])
        if(vt['gen_passlast_eta'][it_p] > -100):
            h_sim_passlast_eta.Fill(vt['gen_eta'][it_p])
           # h_sim_passlast_pt.Fill(vt['gen_pt'][it_p])

        if abs(vt['gen_eta'][it_p]) < 0.9:
            h_sim_pt_barrel.Fill(vt['gen_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 1.4:
            h_sim_pt_trans.Fill(vt['gen_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 2.5:
            h_sim_pt_endcap.Fill(vt['gen_pt'][it_p])

    for it_p in range( vt['np_gen_toReco'][0]): # loop over matched simulated tracks
        h_sim_to_reco_match_eta.Fill(vt['gen_matched_eta'][it_p])

        if abs(vt['gen_eta'][it_p]) < 0.9:
            h_sim_to_reco_match_pt_barrel.Fill(vt['gen_matched_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 1.4:
            h_sim_to_reco_match_pt_trans.Fill(vt['gen_matched_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 2.5:
            h_sim_to_reco_match_pt_endcap.Fill(vt['gen_matched_pt'][it_p])

#efficiency and fake rate wrt eta
h_fakerate_eta = fill_hist_ratio(h_fake_eta, h_reco_eta,"fake_rate_eta")
h_eff_eta = fill_hist_ratio(h_sim_to_reco_match_eta, h_sim_eta, "eff_eta")    

h_eff_hit3_eta = fill_hist_ratio(h_sim_passhit3_eta, h_sim_eta, "eff_hit3_eta")
h_eff_hitlast_eta = fill_hist_ratio(h_sim_passlast_eta, h_sim_eta, "eff_hitlast_eta")

# efficiency and fake rate wrt pt
h_eff_pt_barrel = fill_hist_ratio(h_sim_to_reco_match_pt_barrel, h_sim_pt_barrel, "eff_pt_barrel", binning="log")
h_eff_pt_trans = fill_hist_ratio(h_sim_to_reco_match_pt_trans, h_sim_pt_trans, "eff_pt_trans", binning="log")
h_eff_pt_endcap = fill_hist_ratio(h_sim_to_reco_match_pt_endcap, h_sim_pt_endcap, "eff_pt_endcap", binning="log")

h_fakerate_pt_barrel = fill_hist_ratio(h_fake_pt_barrel, h_reco_pt_barrel, "fake_rate_pt_barrel", binning = "log")
h_fakerate_pt_trans = fill_hist_ratio(h_fake_pt_trans, h_reco_pt_trans, "fake_rate_pt_trans", binning = "log")
h_fakerate_pt_endcap = fill_hist_ratio(h_fake_pt_endcap, h_reco_pt_endcap, "fake_rate_pt_endcap", binning = "log")

outfile = "histograms/trackValHistograms" + args.ptmode  + ".root"    
o = ROOT.TFile(outfile,"recreate")
h_reco_eta.Write()
h_fake_eta.Write()
h_sim_eta.Write()
h_sim_to_reco_match_eta.Write()

h_reco_pt_barrel.Write()
h_reco_pt_endcap.Write()
h_reco_pt_trans.Write()

h_fake_pt_barrel.Write()
h_fake_pt_endcap.Write()
h_fake_pt_trans.Write()

h_sim_pt_barrel.Write()
h_sim_pt_endcap.Write()
h_sim_pt_trans.Write()

h_sim_to_reco_match_pt_barrel.Write()
h_sim_to_reco_match_pt_trans.Write()
h_sim_to_reco_match_pt_endcap.Write()

#-----efficiencies-----------
h_fakerate_eta.Write()
h_eff_eta.Write()
h_eff_hit3_eta.Write()
h_eff_hitlast_eta.Write()

h_fakerate_pt_barrel.Write()
h_fakerate_pt_endcap.Write()
h_fakerate_pt_trans.Write()

h_eff_pt_barrel.Write()
h_eff_pt_endcap.Write()
h_eff_pt_trans.Write()
o.Close()
