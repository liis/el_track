import ROOT, sys
from array import array
from tree_variables import eff_list, var_type, maxhit
from histlib import fill_hist_ratio, fill_hist_ratio_poisson, log_binning, fill_hists_by_eta_regions, initialize_histograms

debug = False
var_list = eff_list

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--testRun', dest='is_test_run', action="store_true", default=False, required=False)
parser.add_argument('--infile', dest='infile', required=True) #file name of the input tree
parser.add_argument('--outdir', dest='outdir', required=False, default="histograms")
args = parser.parse_args()

#indir = "$WORKING_DIR/tree_to_histo/"
indir = ""
infile = args.infile

print "Opening input file: " + indir+infile

f = ROOT.TFile(indir + infile)
t = f.Get("trackValTreeMaker/trackValTree")
#t = f.Get("TrackValTreeMaker/trackValTree") # for old trees (??)

report_every = 10000
if args.is_test_run:
    max_event=1000
else:
    max_event = -1

nEvt = t.GetEntries()

print "Analyzing dataset with " + str(nEvt) + " events"

#### initialize tree #######
vt = dict([ (v, var_type(v)) for v in var_list ]) #associate proper data type for variables in the tree


t.SetBranchStatus("*",0)
print "Performing tree setup..."
for v in var_list: # relate vList to the tree
#    print "starting setting up: " + v
#     print str(vt[v])
    t.SetBranchStatus(v, 1)
    t.SetBranchAddress(v, vt[v])
    t.AddBranchToCache(v, ROOT.kTRUE)

t.StopCacheLearningPhase()

neta = 50
mineta = -2.5
maxeta = 2.5

npt = 50
minpt = 0.1
maxpt = 200

nsimhits = 25
minsimhit = 0
maxsimhit = 25

xbinspt = log_binning(npt,minpt,maxpt) #log binning for pt histograms

eta_regions = ["barrel", "trans", "endcap"]

#seed_vars = ["gen_matchedSeedOkCharge", "gen_matchedSeedQuality"]

simeta_vars = ["sim_eta", "sim_eta_matchedTrack", "sim_eta_matchedSeed", "sim_eta_matchedSeed_matchedCharge"]
simpt_vars = ["sim_pt", "sim_pt_matchedTrack", "sim_pt_matchedSeed", "sim_pt_matchedSeed_matchedCharge"]

recpt_vars = ["reco_pt", "fake_pt"]
receta_vars = ["reco_eta", "fake_eta"]


simeta_hists = initialize_histograms( simeta_vars, bin_reg = (neta, mineta, maxeta) )
simpt_hists = initialize_histograms( simpt_vars, bin_reg = (npt, array('d', xbinspt)), hist_in_regions = eta_regions )

receta_hists = initialize_histograms( receta_vars, bin_reg = (neta, mineta, maxeta) )
recpt_hists = initialize_histograms( recpt_vars, bin_reg = (npt, array('d', xbinspt)), hist_in_regions = eta_regions )

#all_hists = [simeta_hists, simpt_hists, receta_hists, recpt_hists]
all_hists = {
    "sim_eta": simeta_hists,
    "sim_pt": simpt_hists,
    "rec_eta": receta_hists,
    "rec_pt": recpt_hists
    }

# Event loop
for i in range(nEvt):
    if i % report_every == 0: print "Event nr: " + str(i)
    if i == max_event and not i == -1: break

    t.LoadTree(i)
    t.GetEntry(i)
    for it_p in range( vt['np_reco'][0]): # loop over reco-tracks

        reco_eta = vt['reco_eta'][it_p]
        reco_pt = vt['reco_pt'][it_p]

        receta_hists["reco_eta"].Fill(reco_eta)
        if reco_pt > 30:
            fill_hists_by_eta_regions(reco_eta, reco_pt, "reco_pt", recpt_hists)

    for it_p in range( vt['np_fake'][0]): # loop over fake traks
        fake_eta = vt['fake_eta'][it_p]
        fake_pt = vt['fake_pt'][it_p]

        receta_hists["fake_eta"].Fill(fake_eta)
        if fake_pt > 30:
            fill_hists_by_eta_regions(fake_eta, fake_pt, "fake_pt", recpt_hists)

    for it_p in range( vt['np_gen'][0]): #loop over simulated tracks
        gen_track_pt = vt['gen_pt'][it_p]
        gen_track_eta = vt['gen_eta'][it_p]

        simeta_hists["sim_eta"].Fill(gen_track_eta)
        fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt", simpt_hists) # sim pt histograms in 3 eta regions

        # ----------------------------- matched to track ------------------------------
        if vt["gen_reco_matched"][it_p]: # if matched to reco tracks (associatorByHits (matched/rec > 0.75))
            simeta_hists["sim_eta_matchedTrack"].Fill(gen_track_eta)
            fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedTrack", simpt_hists)

        # ------------------------------ matched to seed ------------------------------
        if vt["gen_matchedSeedQuality"][it_p] == 1:
            simeta_hists["sim_eta_matchedSeed"].Fill(vt['gen_eta'][it_p])
            fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed", simpt_hists)
            if vt["gen_matchedSeedOkCharge"][it_p] > 0:
                simeta_hists["sim_eta_matchedSeed_matchedCharge"].Fill(gen_track_eta)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed_matchedCharge", simpt_hists)

eta_regions = ["barrel", "trans", "endcap"]

print "Processing total Sim-to-Reco efficiencies:"
efficiency_histograms={}
efficiency_histograms["h_eff_eta"] = fill_hist_ratio_poisson(simeta_hists["sim_eta_matchedTrack"], simeta_hists["sim_eta"], "eff_eta")
efficiency_histograms["h_eff_seed_eta"] = fill_hist_ratio_poisson(simeta_hists["sim_eta_matchedSeed"], simeta_hists["sim_eta"], "eff_seed_eta")
efficiency_histograms["h_eff_wrt_seed_eta"] = fill_hist_ratio_poisson(simeta_hists["sim_eta_matchedTrack"], simeta_hists["sim_eta_matchedSeed"], "eff_wrt_seed_eta")

for region in eta_regions:
    efficiency_histograms["h_eff_pt" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_" + region], simpt_hists["sim_pt_" + region], "eff_pt_" + region, binning="log")
    efficiency_histograms["h_eff_seed_pt" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedSeed_"+ region], simpt_hists["sim_pt_" + region], "eff_seed_pt_" + region, binning="log")
    efficiency_histograms["h_eff_wrt_seed_pt" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_" + region], simpt_hists["sim_pt_matchedSeed_" + region], "eff_wrt_seed_pt_" + region, binning = "log")


# -------------------- fake rates ---------------------------------
efficiency_histograms["h_fakerate_eta"] = fill_hist_ratio(receta_hists["fake_eta"], receta_hists["reco_eta"],"fake_rate_eta")

h_fakerate_pt = {}
for region in eta_regions:
    efficiency_histograms["h_fakerate_pt_" + region] = fill_hist_ratio(recpt_hists["fake_pt_" + region], recpt_hists["reco_pt_" + region], "fake_rate_pt_" + region, binning = "log")

#---------------------create outfile------------------------
outfilename = (args.infile).split('/trackValTree_',1)[1]
outfile = args.outdir + "/efficiencyHistograms_" + outfilename
print "Writing histograms to file: " + outfile
o = ROOT.TFile(outfile,"recreate")

#-----------------write histograms to file----------------
print "Saving histograms for variables..."
for histos in all_hists: # loop over histogram dictionaries
    dir = o.mkdir(histos)
    dir.cd()
    for hist in all_hists[histos]: # loop over individual histograms
        all_hists[histos][hist].Write()

#------------ write efficiencies and fake rates-----------
print "Saving efficiency histograms..."
dir = o.mkdir("efficiencies")
dir.cd()
for histogram in efficiency_histograms.values():
    histogram.Write()

o.Close()
