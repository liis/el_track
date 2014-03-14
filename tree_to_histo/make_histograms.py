import ROOT, sys
from array import array
from tree_variables import var_list, var_type, maxhit
from histlib import fill_hist_ratio, log_binning, fill_hists_by_eta_regions, initialize_histograms

debug = False

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--testRun', dest='is_test_run', action="store_true", default=False, required=False)
parser.add_argument('--infile', dest='infile', required=True) #file name of the input tree
args = parser.parse_args()

indir = "$WORKING_DIR/tree_to_histo/"

infile = args.infile
print "Opening input file: " + infile

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 1000
if args.is_test_run:
    max_event=1000
else:
    max_event = -1
nEvt = t.GetEntries()

print "Analyzing dataset with " + str(nEvt) + " events"

#### initialize tree #######
vt = dict([ (v, var_type(v)) for v in var_list ]) #associate proper data type for variables in the tree

t.SetBranchStatus("*",0)
for v in var_list: # relate vList to the tree
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

#------------------- initialize histograms--------------------------------------------
#other_vars = ["gen_matched_qoverp", "gen_matched_rec_qoverp", "gen_matched_cotth", "gen_matched_rec_cotth", "gen_matched_phi", "gen_matched_rec_phi", "gen_matched_d0", "gen_matched_rec_d0", "gen_matched_z0 ", "gen_matched_rec_z0" ]

pull_vars = ["pt_pull", "qoverp_pull", "theta_pull", "phi_pull","d0_pull", "z0_pull"]
res_vars = ["pt_res", "qoverp_res", "eta_res", "cotth_res", "phi_res", "d0_res", "z0_res"]

seed_vars = ["gen_matchedSeedOkCharge"]

simeta_vars = ["sim_eta", "sim_eta_smallBrem", "sim_eta_matchedTrack", "sim_eta_matchedTrack_smallBrem", 
               "sim_eta_matchedTrack_tracker", "sim_eta_matchedTrack_ecal", 
               "sim_eta_matchedSeed", "sim_eta_matchedSeed_smallBrem", "sim_eta_matchedSeed_trackerOnly", 
               "sim_eta_matchedSeed_ecalOnly", "sim_eta_matchedSeed_tracker", "sim_eta_matchedSeed_ecal", "sim_eta_matchedSeed_matchedCharge"]
simhit_vars = ["sim_nrhits", "sim_nrhits_matchedSeed", "sim_nrhits_matchedTrack", "sim_nrhits_smallBrem", 
               "sim_nrhits_matchedTrack_smallBrem", "sim_nrhits_matchedSeed_smallBrem", "sim_nrhits_matchedSeed_trackerOnly", 
               "sim_nrhits_matchedSeed_trackerOnly_smallBrem", "sim_nrhits_matchedSeed_ecalOnly", "sim_nrhits_matchedSeed_ecalOnly_smallBrem", 
               "sim_nrhits_matchedSeed_ecal", "sim_nrhits_matchedSeed_tracker"]
simpt_vars = ["sim_pt", "sim_pt_smallBrem", "sim_pt_matchedTrack", "sim_pt_matchedTrack_tracker", "sim_pt_matchedTrack_ecal", 
              "sim_pt_matchedTrack_smallBrem", "sim_pt_matchedSeed", "sim_pt_matchedSeed_tracker", "sim_pt_matchedSeed_trackerOnly",
              "sim_pt_matchedSeed_ecalOnly", "sim_pt_matchedSeed_ecal","sim_pt_matchedSeed_smallBrem", "sim_pt_matchedSeed_matchedCharge"]
recpt_vars = ["reco_pt", "fake_pt"]
receta_vars = ["reco_eta", "fake_eta"]

pull_hists = initialize_histograms( pull_vars, bin_reg = (50, -2, 2))
res_hists = initialize_histograms( res_vars, bin_reg = (100, -1.3, 1.3))

simeta_hists = initialize_histograms( simeta_vars, bin_reg = (neta, mineta, maxeta) )
simhit_hists = initialize_histograms( simhit_vars, bin_reg = (nsimhits, minsimhit, maxsimhit)) 
simhit_byetaregion_hists = initialize_histograms( simhit_vars, bin_reg = (nsimhits, minsimhit, maxsimhit), hist_in_regions = eta_regions)
simpt_hists = initialize_histograms( simpt_vars, bin_reg = (npt, array('d', xbinspt)), hist_in_regions = eta_regions )

receta_hists = initialize_histograms( receta_vars, bin_reg = (neta, mineta, maxeta) )
recpt_hists = initialize_histograms( recpt_vars, bin_reg = (npt, array('d', xbinspt)), hist_in_regions = eta_regions )

all_hists = [simeta_hists, simhit_hists, simhit_byetaregion_hists, simpt_hists, receta_hists, recpt_hists, pull_hists, res_hists]

#------------------------initialize tracker-hit-histograms---------------------
sim_hit_quality_flags = {"quality065": 0.65,"quality075": 0.75,"quality085": 0.85, "quality095": 0.95}

h_sim_eta_byhit = []
h_sim_pt_byhit_barrel = []
h_sim_pt_byhit_trans = []
h_sim_pt_byhit_endcap = []

h_sim_eta_byhit_quality = {}
h_sim_pt_byhit_quality_barrel = {}
h_sim_pt_byhit_quality_trans = {}
h_sim_pt_byhit_quality_endcap = {}


for nrhit in range(0,maxhit):
    h_sim_eta_byhit.append(ROOT.TH1F("eta_at_"+str(nrhit), "eta_at_"+str(nrhit),neta,mineta,maxeta) )
    
    h_sim_pt_byhit_barrel.append(ROOT.TH1F("pt_at_"+str(nrhit)+"_barrel", "pt_at_"+str(nrhit),npt,array('d',xbinspt)) )
    h_sim_pt_byhit_trans.append(ROOT.TH1F("pt_at_"+str(nrhit)+"_trans", "pt_at_"+str(nrhit),npt, array('d', xbinspt)) )
    h_sim_pt_byhit_endcap.append(ROOT.TH1F("pt_at_"+str(nrhit)+"_endcap", "pt_at_"+str(nrhit),npt, array('d', xbinspt)) )


for quality_flag in sim_hit_quality_flags:
    h_sim_eta_byhit_quality[quality_flag] = []
    
    h_sim_pt_byhit_quality_barrel[quality_flag] = []
    h_sim_pt_byhit_quality_endcap[quality_flag] = []
    h_sim_pt_byhit_quality_trans[quality_flag] = []
    for nrhit in range(0,maxhit):
        h_sim_eta_byhit_quality[quality_flag].append(ROOT.TH1F("eta_at_"+str(nrhit)+"_"+quality_flag,"eta_at_"+str(nrhit)+"_"+quality_flag, neta, mineta, maxeta) ) #pt for tracks, that have nth track (save for efficiencies)
        h_sim_pt_byhit_quality_barrel[quality_flag].append(ROOT.TH1F("pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", "pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", npt, array('d', xbinspt)) )
        h_sim_pt_byhit_quality_endcap[quality_flag].append(ROOT.TH1F("pt_at"+str(nrhit)+"_"+quality_flag+"_endcap", "pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", npt, array('d', xbinspt)) )
        h_sim_pt_byhit_quality_trans[quality_flag].append(ROOT.TH1F("pt_at"+str(nrhit)+"_"+quality_flag+"_trans", "pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", npt, array('d',xbinspt)) )
#--------------------------------------------------------------------


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
        fill_hists_by_eta_regions(reco_eta, reco_pt, "reco_pt", recpt_hists)

    for it_p in range( vt['np_fake'][0]): # loop over fake traks
        fake_eta = vt['reco_eta'][it_p]
        fake_pt = vt['fake_eta'][it_p]

        receta_hists["fake_eta"].Fill(fake_eta)
        fill_hists_by_eta_regions(fake_eta, fake_pt, "fake_pt", recpt_hists)

    for it_p in range( vt['np_gen'][0]): #loop over simulated tracks

        gen_track_nrSimHits = vt['gen_nrUniqueSimHits'][it_p]        
        gen_track_pt = vt['gen_pt'][it_p]
        gen_track_eta = vt['gen_eta'][it_p]

        simeta_hists["sim_eta"].Fill(gen_track_eta)

        fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt", simpt_hists) # sim pt histograms in 3 eta regions

        simhit_hists["sim_nrhits"].Fill(gen_track_nrSimHits)        
        fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits", simhit_byetaregion_hists)

        if vt['gen_bremFraction'][it_p] < 0.2: # with reasonable brem
            simeta_hists["sim_eta_smallBrem"].Fill(gen_track_eta)
            fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_smallBrem", simpt_hists)
            
            simhit_hists["sim_nrhits_smallBrem"].Fill(gen_track_nrSimHits)
            fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_smallBrem", simhit_byetaregion_hists)

        if vt["gen_reco_matched"][it_p]: # if matched to reco tracks
            #---------------pull ja resolution histos------------
            pull_hists["pt_pull"].Fill(vt["pt_pull"][it_p])
            pull_hists["qoverp_pull"].Fill(vt["qoverp_pull"][it_p])
            pull_hists["theta_pull"].Fill(vt["theta_pull"][it_p])
            pull_hists["phi_pull"].Fill(vt["phi_pull"][it_p])
            pull_hists["d0_pull"].Fill(vt["d0_pull"][it_p])
            pull_hists["z0_pull"].Fill(vt["z0_pull"][it_p])
            
#            res_hists["eta_res"].Fill(vt["gen_matched_eta"][it_p] - vt["gen_matched_eta"][it_p]) #next run
            res_hists["pt_res"].Fill( (vt["gen_matched_rec_pt"][it_p] - vt["gen_matched_pt"][it_p])/vt["gen_matched_pt"][it_p] )
            res_hists["qoverp_res"].Fill( vt["gen_matched_rec_qoverp"][it_p] - vt["gen_matched_qoverp"][it_p] )
            res_hists["phi_res"].Fill(vt["gen_matched_rec_phi"][it_p] - vt["gen_matched_phi"][it_p])
            res_hists["cotth_res"].Fill(vt["gen_matched_rec_cotth"][it_p] - vt["gen_matched_cotth"][it_p])
            res_hists["z0_res"].Fill(vt["gen_matched_rec_z0"][it_p] - vt["gen_matched_z0"][it_p])
            res_hists["d0_res"].Fill(vt["gen_matched_rec_d0"][it_p] - vt["gen_matched_d0"][it_p])
                                        
            #---------------------------------------------------
            simhit_hists["sim_nrhits_matchedTrack"].Fill(gen_track_nrSimHits)        
            fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedTrack", simhit_byetaregion_hists)

            simeta_hists["sim_eta_matchedTrack"].Fill(gen_track_eta)
            fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedTrack", simpt_hists)
            if vt["is_ecalDrivenSeed"][it_p] == 1:
                simeta_hists["sim_eta_matchedTrack_ecal"].Fill(gen_track_eta)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedTrack_ecal", simpt_hists)
            if vt["is_trackerDrivenSeed"][it_p] == 1:
                simeta_hists["sim_eta_matchedTrack_tracker"].Fill(gen_track_eta)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedTrack_tracker", simpt_hists)

            if vt['gen_bremFraction'][it_p] < 0.2: 
                simhit_hists["sim_nrhits_matchedTrack_smallBrem"].Fill(gen_track_nrSimHits)
                simeta_hists["sim_eta_matchedTrack_smallBrem"].Fill(gen_track_eta)

                fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedTrack_smallBrem", simpt_hists)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedTrack_smallBrem", simhit_byetaregion_hists)

#        print len(vt["gen_matchedSeedOkCharge"])      
#        print "seed charge match = " + str(vt["gen_matchedSeedOkCharge"][it_p]) 
        if vt['gen_nrMatchedSeedHits'][it_p] > 1: # if matched to seed
            simeta_hists["sim_eta_matchedSeed"].Fill(vt['gen_eta'][it_p])
            fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed", simpt_hists)
            
            try:
                if vt["gen_matchedSeedOkCharge"][it_p] > 0:
                    simeta_hists["sim_eta_matchedSeed_matchedCharge"].Fill(gen_track_eta)
                    fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed_matchedCharge", simpt_hists)
            except IndexError:
                print "Skip filling seed charge efficiency --> FIX BUG!"

            simhit_hists["sim_nrhits_matchedSeed"].Fill(vt['gen_nrUniqueSimHits'][it_p])
            fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedSeed", simhit_byetaregion_hists)
            if vt['gen_bremFraction'][it_p] < 0.2:
                simeta_hists["sim_eta_matchedSeed_smallBrem"].Fill(gen_track_eta)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed_smallBrem", simpt_hists)

                simhit_hists["sim_nrhits_matchedSeed_smallBrem"].Fill(gen_track_nrSimHits)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedSeed_smallBrem", simhit_byetaregion_hists)
                if vt['is_ecalDrivenSeed'][it_p] == 1 and vt["is_trackerDrivenSeed"][it_p]== 0:
                    simhit_hists["sim_nrhits_matchedSeed_ecalOnly_smallBrem"].Fill(gen_track_nrSimHits)
                    fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedSeed_ecalOnly_smallBrem", simhit_byetaregion_hists)
                elif vt['is_ecalDrivenSeed'][it_p] == 0 and vt["is_trackerDrivenSeed"][it_p] == 1:
                    simhit_hists["sim_nrhits_matchedSeed_trackerOnly_smallBrem"].Fill(gen_track_nrSimHits)
                    fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedSeed_trackerOnly_smallBrem", simhit_byetaregion_hists)

            if vt['is_ecalDrivenSeed'][it_p] == 1:
                simhit_hists["sim_nrhits_matchedSeed_ecal"].Fill(gen_track_nrSimHits)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedSeed_ecal", simhit_byetaregion_hists)
                
                simeta_hists["sim_eta_matchedSeed_ecal"].Fill(gen_track_eta)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed_ecal", simpt_hists)

                if vt["is_trackerDrivenSeed"][it_p] == 0:
                    simhit_hists["sim_nrhits_matchedSeed_ecalOnly"].Fill(gen_track_nrSimHits)
                    simeta_hists["sim_eta_matchedSeed_ecalOnly"].Fill(gen_track_eta)
                    fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed_ecalOnly", simpt_hists)
                    fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedSeed_ecalOnly", simhit_byetaregion_hists) 
            if vt['is_trackerDrivenSeed'][it_p] == 1:
                simhit_hists["sim_nrhits_matchedSeed_tracker"].Fill(gen_track_nrSimHits)
                fill_hists_by_eta_regions(gen_track_eta,gen_track_nrSimHits, "sim_nrhits_matchedSeed_tracker", simhit_byetaregion_hists)
                simeta_hists["sim_eta_matchedSeed_tracker"].Fill(gen_track_eta)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed_tracker", simpt_hists)

                if vt["is_ecalDrivenSeed"][it_p] == 0:
                    simhit_hists["sim_nrhits_matchedSeed_trackerOnly"].Fill(gen_track_nrSimHits)
                    simeta_hists["sim_eta_matchedSeed_trackerOnly"].Fill(gen_track_eta)
                    fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt_matchedSeed_trackerOnly", simpt_hists)
                    fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits_matchedSeed_trackerOnly", simhit_byetaregion_hists)
################### Analyze sim hits ##############################
        for nrhit in range(0,maxhit):
            pt_at_entry = vt["gen_hit_pt"][it_p][nrhit]
            hit_subdet = vt["gen_hit_subdetector"][it_p][nrhit]
            hit_layer = vt["gen_hit_layer"][it_p][nrhit]

#            print "hit pt = " + str(pt_at_entry) + ", hit subdetector = " + str(hit_subdet) + ", layer = " + str(hit_layer);
            if pt_at_entry > 0: # Nth hit exists -- fill eta and pt histograms for efficiency plots
                h_sim_eta_byhit[nrhit].Fill(gen_track_eta)
                
                if abs(vt['gen_eta'][it_p]) < 0.9:
                    h_sim_pt_byhit_barrel[nrhit].Fill(gen_track_pt)
                elif abs(vt['gen_eta'][it_p]) < 1.4:
                    h_sim_pt_byhit_trans[nrhit].Fill(gen_track_pt)
                elif abs(vt['gen_eta'][it_p]) < 2.5:
                    h_sim_pt_byhit_endcap[nrhit].Fill(gen_track_pt)

                for quality_flag, cut in sim_hit_quality_flags.iteritems():
                    if pt_at_entry/gen_track_pt > cut: # check hit pt-quality condition at each hit 
                        h_sim_eta_byhit_quality[quality_flag][nrhit].Fill(gen_track_eta)

                        if abs(vt['gen_eta'][it_p]) < 0.9:
                            h_sim_pt_byhit_quality_barrel[quality_flag][nrhit].Fill(gen_track_pt) 
                        elif abs(vt['gen_eta'][it_p]) < 1.4:
                            h_sim_pt_byhit_quality_trans[quality_flag][nrhit].Fill(gen_track_pt) 
                        elif abs(vt['gen_eta'][it_p]) < 2.5:
                            h_sim_pt_byhit_quality_endcap[quality_flag][nrhit].Fill(gen_track_pt) 

#----------------------------------hit-by-hit efficiencies---------------------------------
h_eff_hit_eta = {}
h_eff_hit_pt_barrel = {}
h_eff_hit_pt_endcap = {}
h_eff_hit_pt_trans = {}
for quality_flag in sim_hit_quality_flags:
    h_eff_hit_eta[quality_flag] = []
    h_eff_hit_pt_barrel[quality_flag] = []
    h_eff_hit_pt_trans[quality_flag] = []
    h_eff_hit_pt_endcap[quality_flag] = []

    for nrhit in range(0,maxhit):
        h_eff_hit_eta[quality_flag].append( fill_hist_ratio(h_sim_eta_byhit_quality[quality_flag][nrhit], h_sim_eta_byhit[nrhit], "eff_eta_athit_"+str(nrhit)+"_"+quality_flag) )
        
        h_eff_hit_pt_barrel[quality_flag].append( fill_hist_ratio(h_sim_pt_byhit_quality_barrel[quality_flag][nrhit], h_sim_pt_byhit_barrel[nrhit], "eff_pt_athit_"+str(nrhit)+"_barrel_"+quality_flag, binning = "log") )
        h_eff_hit_pt_trans[quality_flag].append( fill_hist_ratio(h_sim_pt_byhit_quality_trans[quality_flag][nrhit], h_sim_pt_byhit_trans[nrhit], "eff_pt_athit_"+str(nrhit)+"_trans_"+quality_flag, binning = "log") )
        h_eff_hit_pt_endcap[quality_flag].append( fill_hist_ratio(h_sim_pt_byhit_quality_endcap[quality_flag][nrhit], h_sim_pt_byhit_endcap[nrhit], "eff_pt_athit_"+str(nrhit)+"_endcap_"+quality_flag, binning = "log") )

nbin = 15
h_nrhit_barrel = {}
h_nrhit_endcap = {}
h_nrhit_trans = {}

for quality_flag in sim_hit_quality_flags:
    h_nrhit_barrel[quality_flag] = ROOT.TH1F("nr_hit_" + quality_flag + "_barrel", "nr_hit_" + quality_flag + "_barrel", nbin, 1, nbin+1)
    h_nrhit_trans[quality_flag] = ROOT.TH1F("nr_hit_" + quality_flag + "_trans", "nr_hit_" + quality_flag + "_trans", nbin, 1, nbin+1)
    h_nrhit_endcap[quality_flag] = ROOT.TH1F("nr_hit_" + quality_flag + "_endcap", "nr_hit_" + quality_flag + "_endcap", nbin, 1, nbin+1)

    for ibin in range(1,nbin):
        if h_sim_pt_byhit_barrel[ibin-1].Integral(): h_nrhit_barrel[quality_flag].SetBinContent(ibin, h_sim_pt_byhit_quality_barrel[quality_flag][ibin-1].Integral()/h_sim_pt_byhit_barrel[ibin-1].Integral())
        if h_sim_pt_byhit_trans[ibin-1].Integral(): h_nrhit_trans[quality_flag].SetBinContent(ibin, h_sim_pt_byhit_quality_trans[quality_flag][ibin-1].Integral()/h_sim_pt_byhit_trans[ibin-1].Integral())
        if h_sim_pt_byhit_endcap[ibin-1].Integral(): h_nrhit_endcap[quality_flag].SetBinContent(ibin, h_sim_pt_byhit_quality_endcap[quality_flag][ibin-1].Integral()/h_sim_pt_byhit_endcap[ibin-1].Integral())
    

eta_regions = ["barrel", "trans", "endcap"]

print "Processing total Sim-to-Reco efficiencies:"
efficiency_histograms={}

efficiency_histograms["h_eff_eta"] = fill_hist_ratio(simeta_hists["sim_eta_matchedTrack"], simeta_hists["sim_eta"], "eff_eta")
efficiency_histograms["h_eff_eta_smallBrem"] = fill_hist_ratio(simeta_hists["sim_eta_matchedTrack_smallBrem"], simeta_hists["sim_eta_smallBrem"], "eff_eta_smallBrem")

efficiency_histograms["h_eff_seed_eta"] = fill_hist_ratio(simeta_hists["sim_eta_matchedSeed"], simeta_hists["sim_eta"], "eff_seed_eta")
efficiency_histograms["h_eff_seed_eta_smallBrem"] = fill_hist_ratio(simeta_hists["sim_eta_matchedSeed_smallBrem"], simeta_hists["sim_eta_smallBrem"], "eff_seed_eta_smallBrem")
efficiency_histograms["h_eff_seed_eta_tracker"] = fill_hist_ratio(simeta_hists["sim_eta_matchedSeed_tracker"], simeta_hists["sim_eta"], "eff_seed_eta_tracker")
efficiency_histograms["h_eff_seed_eta_ecal"] =  fill_hist_ratio(simeta_hists["sim_eta_matchedSeed_ecal"], simeta_hists["sim_eta"], "eff_seed_eta_ecal")

efficiency_histograms["h_eff_wrt_seed_eta"] = fill_hist_ratio(simeta_hists["sim_eta_matchedTrack"], simeta_hists["sim_eta_matchedSeed"], "eff_wrt_seed_eta")
efficiency_histograms["h_eff_wrt_seed_eta_smallBrem"] = fill_hist_ratio(simeta_hists["sim_eta_matchedTrack_smallBrem"], simeta_hists["sim_eta_matchedSeed_smallBrem"], "eff_wrt_seed_eta_smallBrem")
efficiency_histograms["h_eff_wrt_seed_eta_tracker"] = fill_hist_ratio(simeta_hists["sim_eta_matchedTrack_tracker"], simeta_hists["sim_eta_matchedSeed_tracker"], "eff_wrt_seed_eta_tracker")
efficiency_histograms["h_eff_wrt_seed_eta_ecal"] = fill_hist_ratio(simeta_hists["sim_eta_matchedTrack_ecal"], simeta_hists["sim_eta_matchedSeed"], "eff_wrt_seed_eta_ecal")

efficiency_histograms["h_eff_nrhits"] = fill_hist_ratio(simhit_hists["sim_nrhits_matchedTrack"], simhit_hists["sim_nrhits"], "eff_nrhits")
efficiency_histograms["h_eff_wrt_seed_nrhits"] = fill_hist_ratio(simhit_hists["sim_nrhits_matchedTrack"], simhit_hists["sim_nrhits_matchedSeed"], "eff_wrt_seed_nrhits")
efficiency_histograms["h_eff_seed_nrhits"] = fill_hist_ratio(simhit_hists["sim_nrhits_matchedSeed"], simhit_hists["sim_nrhits"], "eff_seed_nrhits")

efficiency_histograms["h_eff_nrhits_smallBrem"] = fill_hist_ratio(simhit_hists["sim_nrhits_matchedTrack_smallBrem"], simhit_hists["sim_nrhits_smallBrem"], "eff_nrhits_smallBrem")
efficiency_histograms["h_eff_wrt_seed_nrhits_smallBrem"] = fill_hist_ratio(simhit_hists["sim_nrhits_matchedTrack_smallBrem"], simhit_hists["sim_nrhits_matchedSeed_smallBrem"], "eff_wrt_seed_nrhits_smallBrem")
efficiency_histograms["h_eff_seed_nrhits_smallBrem"] = fill_hist_ratio(simhit_hists["sim_nrhits_matchedSeed_smallBrem"], simhit_hists["sim_nrhits_smallBrem"], "eff_seed_nrhits_smallBrem")
efficiency_histograms["h_eff_seed_charge_eta"] = fill_hist_ratio(simeta_hists["sim_eta_matchedSeed_matchedCharge"], simeta_hists["sim_eta_matchedSeed"], "eff_seed_charge_eta")

for region in eta_regions:
    efficiency_histograms["h_eff_pt" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_" + region], simpt_hists["sim_pt_" + region], "eff_pt_" + region, binning="log")
    efficiency_histograms["h_eff_pt_smallBrem" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_smallBrem_" + region], simpt_hists["sim_pt_smallBrem_" + region], "eff_pt_smallBrem_" + region, binning="log")

    efficiency_histograms["h_eff_seed_pt" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedSeed_"+ region], simpt_hists["sim_pt_" + region], "eff_seed_pt_" + region, binning="log")
    efficiency_histograms["h_eff_seed_pt_smallBrem"+ region] = fill_hist_ratio(simpt_hists["sim_pt_matchedSeed_smallBrem_" + region], simpt_hists["sim_pt_smallBrem_" + region], "eff_seed_pt_smallBrem_" + region, binning="log")
    efficiency_histograms["h_eff_seed_pt_ecal" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedSeed_ecal_"+ region], simpt_hists["sim_pt_" + region], "eff_seed_pt_ecal_" + region, binning="log")
    efficiency_histograms["h_eff_seed_pt_tracker" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedSeed_tracker_"+ region], simpt_hists["sim_pt_" + region], "eff_seed_pt_tracker_" + region, binning="log")
    efficiency_histograms["h_eff_seed_charge_pt" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedSeed_matchedCharge_" + region], simpt_hists["sim_pt_matchedSeed_" + region], "eff_seed_charge_pt_" + region, binning="log")

    efficiency_histograms["h_eff_wrt_seed_pt_tracker" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_tracker_"+ region], simpt_hists["sim_pt_matchedSeed_tracker_" + region], "eff_wrt_seed_pt_tracker_" + region, binning="log")

    efficiency_histograms["h_eff_wrt_seed_pt_ecal" + region] =fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_ecal_"+ region], simpt_hists["sim_pt_matchedSeed_ecal_" + region], "eff_wrt_seed_pt_ecal_" + region, binning="log")

    efficiency_histograms["h_eff_wrt_seed_pt" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_" + region], simpt_hists["sim_pt_matchedSeed_" + region], "eff_wrt_seed_pt_" + region, binning = "log")
    efficiency_histograms["h_eff_wrt_seed_pt_smallBrem" + region] = fill_hist_ratio(simpt_hists["sim_pt_matchedTrack_smallBrem_" + region], simpt_hists["sim_pt_matchedSeed_smallBrem_" + region], "eff_wrt_seed_pt_smallBrem_" + region, binning="log")

    efficiency_histograms["h_eff_nrhits_byetaregion" + region] = fill_hist_ratio(simhit_byetaregion_hists["sim_nrhits_matchedTrack_" + region], simhit_byetaregion_hists["sim_nrhits_" + region], "eff_nhits_" + region)
    efficiency_histograms["h_eff_seed_nrhits_byetaregion" + region] = fill_hist_ratio(simhit_byetaregion_hists["sim_nrhits_matchedSeed_" + region], simhit_byetaregion_hists["sim_nrhits_" + region], "eff_seed_nhits_" + region)
    efficiency_histograms["h_eff_wrt_seed_nrhits_byetaregion" + region] = fill_hist_ratio(simhit_byetaregion_hists["sim_nrhits_matchedTrack_" + region], simhit_byetaregion_hists["sim_nrhits_matchedSeed_" + region], "eff_wrt_seed_nhits_" + region)
  
    efficiency_histograms["h_eff_nrhits_byetaregion_smallBrem" + region] = fill_hist_ratio(simhit_byetaregion_hists["sim_nrhits_matchedTrack_smallBrem_" + region], simhit_byetaregion_hists["sim_nrhits_smallBrem_" + region], "eff_nhits_smallBrem_" + region)
    efficiency_histograms["h_eff_seed_nrhits_byetaregion_smallBrem" + region] = fill_hist_ratio(simhit_byetaregion_hists["sim_nrhits_matchedSeed_smallBrem_" + region], simhit_byetaregion_hists["sim_nrhits_smallBrem_" + region], "eff_seed_nhits_smallBrem_" + region)
    efficiency_histograms["h_eff_wrt_seed_nrhits_byetaregion_smallBrem" + region] = fill_hist_ratio(simhit_byetaregion_hists["sim_nrhits_matchedTrack_smallBrem_" + region], simhit_byetaregion_hists["sim_nrhits_matchedSeed_smallBrem_" + region], "eff_wrt_seed_nhits_smallBrem_" + region)    

# -------------------- fake rates ---------------------------------
efficiency_histograms["h_fakerate_eta"] = fill_hist_ratio(receta_hists["fake_eta"], receta_hists["reco_eta"],"fake_rate_eta")
h_fakerate_pt = {}
for region in eta_regions:
    efficiency_histograms["h_fakerate_pt_" + region] = fill_hist_ratio(recpt_hists["fake_pt_" + region], recpt_hists["reco_pt_" + region], "fake_rate_pt_" + region, binning = "log")

#---------------------create outfile------------------------
outfilename = (args.infile).split('/trackValTree_',1)[1]
outfile = "histograms/trackValHistograms_" + outfilename

print "Writing histograms to file: " + outfile
o = ROOT.TFile(outfile,"recreate")

#-----------------write histograms to file----------------
for histos in all_hists: # loop over histogram dictionaries
    for hist in histos: # loop over individual histograms
        histos[hist].Write()

#------------ write efficiencies and fake rates-----------
for histogram in efficiency_histograms.values():
    histogram.Write()

#------------hit-by-hit comparison------------
dir = o.mkdir("VariablesBySimhit")
dir.cd()
for quality_flag in sim_hit_quality_flags:
    h_nrhit_barrel[quality_flag].Write()
    h_nrhit_trans[quality_flag].Write()
    h_nrhit_endcap[quality_flag].Write()
    for nrhit in range(0,maxhit):
        h_eff_hit_eta[quality_flag][nrhit].Write()
        h_sim_eta_byhit_quality[quality_flag][nrhit].Write()

        h_eff_hit_pt_barrel[quality_flag][nrhit].Write()
        h_eff_hit_pt_endcap[quality_flag][nrhit].Write()
        h_eff_hit_pt_trans[quality_flag][nrhit].Write()

        h_sim_pt_byhit_quality_barrel[quality_flag][nrhit].Write()
        h_sim_pt_byhit_quality_endcap[quality_flag][nrhit].Write()
        h_sim_pt_byhit_quality_trans[quality_flag][nrhit].Write()
for nrhit in range(0,maxhit):
    h_sim_eta_byhit[nrhit].Write()

    h_sim_pt_byhit_barrel[nrhit].Write()
    h_sim_pt_byhit_endcap[nrhit].Write()
    h_sim_pt_byhit_trans[nrhit].Write()

o.Close()

