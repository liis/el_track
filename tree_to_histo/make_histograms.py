import ROOT, sys
from array import array
from tree_variables import var_list, var_type, maxhit
from histlib import fill_hist_ratio, log_binning, fill_hists_by_eta_regions, initialize_histograms, initialize_hist_byEtaReg

debug = False

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ptmode', dest='ptmode', choices=["Pt1","Pt10","Pt100","FlatPt"], required=True, help= "pt cut in analyzed dataset")
parser.add_argument('--isGSF', dest='is_gsf',  action="store_true")
args = parser.parse_args()

indir = "$WORKING_DIR/tree_to_histo/input_trees/"

if args.is_gsf:
    infile = "trackValTree_el" + args.ptmode + "_GSF.root"
else:
    infile = "trackValTree_el" + args.ptmode + ".root"

print "Opening input file: " + infile

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 100
max_event = 10000
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

npt = 50
minpt = 0.1
maxpt = 200

nsimhits = 25
minsimhit = 0
maxsimhit = 25

xbinspt = log_binning(npt,minpt,maxpt) #log binning for pt histograms  

eta_regions = ["barrel", "trans", "endcap"]

#------------------- initialize histograms--------------------------------------------

simeta_vars = ["sim_eta", "sim_eta_trackerOnly", "sim_eta_ecalOnly", "sim_eta_matchedTrack", "sim_eta_matchedSeed", "sim_eta_matchedSeed_trackerOnly", "sim_eta_matchedSeed_ecalOnly"]
simhit_vars = ["sim_nrhits", "sim_nrhits_matchedSeed", "sim_nrhits_matchedTrack", "sim_nrhits_smallBrem", "sim_nrhits_matchedTrack_smallBrem", "sim_nrhits_matchedSeed_smallBrem", "sim_nrhits_matchedSeed_trackerOnly_smallBrem", "sim_nrhits_matchedSeed_ecalOnly", "sim_nrhits_matchedSeed_ecalOnly_smallBrem"]
simpt_vars = ["sim_pt", "sim_pt_matchedTrack", "sim_pt_matchedSeed"]

recpt_vars = ["reco_pt", "fake_pt"]
receta_vars = ["reco_eta", "fake_eta"]

simeta_hists = initialize_histograms( simeta_vars, bin_reg = (neta, mineta, maxeta) )
simhit_hists = initialize_histograms( simhit_vars, bin_reg = (nsimhits, minsimhit, maxsimhit)) 
simhit_byetaregion_hists = initialize_histograms( simhit_vars, bin_reg = (nsimhits, minsimhit, maxsimhit), hist_in_regions = eta_regions)
simpt_hists = initialize_histograms( simpt_vars, bin_reg = (npt, array('d', xbinspt)), hist_in_regions = eta_regions )

receta_hists = initialize_histograms( receta_vars, bin_reg = (neta, mineta, maxeta) )
recpt_hists = initialize_histograms( recpt_vars, bin_reg = (npt, array('d', xbinspt)), hist_in_regions = eta_regions )

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
        if vt['is_ecalDrivenSeed'][it_p] == 1 and vt["is_trackerDrivenSeed"][it_p]== 0:
            simeta_hists["sim_eta_ecalOnly"].Fill(gen_track_eta)
        elif vt['is_ecalDrivenSeed'][it_p] == 0 and vt["is_trackerDrivenSeed"][it_p] == 1:
            simeta_hists["sim_eta_trackerOnly"].Fill(gen_track_eta)

        fill_hists_by_eta_regions(gen_track_eta, gen_track_pt, "sim_pt", simpt_hists) # sim pt histograms in 3 eta regions

        simhit_hists["sim_nrhits"].Fill(gen_track_nrSimHits)        
        fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, "sim_nrhits", simhit_byetaregion_hists)

        if vt['gen_bremFraction'][it_p] < 0.2: # with reasonable brem
            h_sim_nrhits_smallBrem.Fill(gen_track_nrSimHits)
            fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, h_sim_nrhits_byetareg_smallBrem)

        if vt["gen_reco_matched"][it_p]: # if matched to reco tracks
            h_sim_nrhits_matchedTrack.Fill(gen_track_nrSimHits)        
            fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, h_sim_nrhits_byetareg_matchedTrack)
            if vt['gen_bremFraction'][it_p] < 0.2: 
                h_sim_nrhits_matchedTrack_smallBrem.Fill(gen_track_nrSimHits)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, h_sim_nrhits_byetareg_matchedTrack_smallBrem)

        if vt['gen_nrMatchedSeedHits'][it_p] > 1: # if matched to seed
            h_sim_eta_matchedSeed.Fill(vt['gen_eta'][it_p])
            fill_hists_by_eta_regions(vt['gen_eta'][it_p], vt['gen_pt'][it_p], h_sim_pt_matchedSeed)

            h_sim_nrhits_matchedSeed.Fill(vt['gen_nrUniqueSimHits'][it_p])
            fill_hists_by_eta_regions(vt['gen_eta'][it_p], vt['gen_nrUniqueSimHits'][it_p], h_sim_nrhits_byetareg_matchedSeed)
            if vt['gen_bremFraction'][it_p] < 0.2:
                h_sim_nrhits_matchedSeed_smallBrem.Fill(gen_track_nrSimHits)
                fill_hists_by_eta_regions(gen_track_eta, gen_track_nrSimHits, h_sim_nrhits_byetareg_matchedSeed_smallBrem)
                if vt['is_ecalDrivenSeed'][it_p] == 1 and vt["is_trackerDrivenSeed"][it_p]== 0:
                    h_sim_nrhits_matchedSeed_smallBrem_ecalOnly.Fill(gen_track_nrSimHits)
                elif vt['is_ecalDrivenSeed'][it_p] == 0 and vt["is_trackerDrivenSeed"][it_p] == 1:
                    h_sim_nrhits_matchedSeed_smallBrem_trackerOnly.Fill(gen_track_nrSimHits)

            if vt['is_ecalDrivenSeed'][it_p] == 1 and vt["is_trackerDrivenSeed"][it_p] == 0:
                h_sim_nrhits_matchedSeed_ecalOnly.Fill(gen_track_nrSimHits)
                h_sim_eta_matchedSeed_ecalOnly.Fill(gen_track_eta)
            elif vt['is_ecalDrivenSeed'][it_p] == 0 and vt["is_trackerDrivenSeed"][it_p] == 1:
                h_sim_nrhits_matchedSeed_trackerOnly.Fill(gen_track_nrSimHits)
                h_sim_eta_matchedSeed_trackerOnly.Fill(gen_track_eta)
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

#############################################################
                
#        if abs(vt['gen_eta'][it_p]) < 0.9:
#            h_sim_pt_barrel.Fill(vt['gen_pt'][it_p])
#        elif abs(vt['gen_eta'][it_p]) < 1.4:
#            h_sim_pt_trans.Fill(vt['gen_pt'][it_p])
#        elif abs(vt['gen_eta'][it_p]) < 2.5:
#            h_sim_pt_endcap.Fill(vt['gen_pt'][it_p])
            

    for it_pm in range( vt['np_gen_toReco'][0]): # loop over matched simulated tracks
        h_sim_to_reco_match_eta.Fill(vt['gen_matched_eta'][it_pm]) 
        fill_hists_by_eta_regions(vt['gen_eta'][it_pm], vt['gen_pt'][it_pm], h_sim_to_reco_match_pt) # sim pt histograms in 3 eta regions


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
        h_nrhit_barrel[quality_flag].SetBinContent(ibin, h_sim_pt_byhit_quality_barrel[quality_flag][ibin-1].Integral()/h_sim_pt_byhit_barrel[ibin-1].Integral())
        h_nrhit_trans[quality_flag].SetBinContent(ibin, h_sim_pt_byhit_quality_trans[quality_flag][ibin-1].Integral()/h_sim_pt_byhit_trans[ibin-1].Integral())
        h_nrhit_endcap[quality_flag].SetBinContent(ibin, h_sim_pt_byhit_quality_endcap[quality_flag][ibin-1].Integral()/h_sim_pt_byhit_endcap[ibin-1].Integral())
    

eta_regions = ["barrel", "trans", "endcap"]
# --------------------------- efficiencies ------------------------
print "Processing total Sim-to-Reco efficiencies:"
h_eff_eta = fill_hist_ratio(h_sim_to_reco_match_eta, h_sim_eta, "eff_eta")
h_eff_seed_eta = fill_hist_ratio(h_sim_eta_matchedSeed, h_sim_eta, "eff_seed_eta")
h_eff_seed_eta_ecalOnly = fill_hist_ratio(h_sim_eta_matchedSeed_ecalOnly, h_sim_eta_ecalOnly, "eff_seed_eta_ecalOnly")
h_eff_seed_eta_trackerOnly = fill_hist_ratio(h_sim_eta_matchedSeed_trackerOnly, h_sim_eta_trackerOnly, "eff_seed_eta_trackerOnly")

h_eff_wrt_seed_eta = fill_hist_ratio(h_sim_to_reco_match_eta, h_sim_eta_matchedSeed, "eff_wrt_seed_eta")

h_eff_nrhits = fill_hist_ratio(h_sim_nrhits_matchedTrack, h_sim_nrhits, "eff_nrhits")
h_eff_seed_nrhits = fill_hist_ratio(h_sim_nrhits_matchedSeed, h_sim_nrhits, "eff_seed_nrhits")

h_eff_nrhits_smallBrem = fill_hist_ratio(h_sim_nrhits_matchedTrack_smallBrem, h_sim_nrhits_smallBrem, "eff_nrhits_smallBrem")
h_eff_seed_nrhits_smallBrem = fill_hist_ratio(h_sim_nrhits_matchedSeed_smallBrem, h_sim_nrhits_smallBrem, "eff_seed_nrhits_smallBrem")

h_eff_pt = {}
h_eff_seed_pt = {}
h_eff_wrt_seed_pt = {}

h_eff_nrhits_byetaregion = {}
h_eff_seed_nrhits_byetaregion = {}

h_eff_nrhits_byetaregion_smallBrem = {}
h_eff_seed_nrhits_byetaregion_smallBrem = {}
for region in eta_regions:
    h_eff_pt[region] = fill_hist_ratio(h_sim_to_reco_match_pt[region], h_sim_pt[region], "eff_pt_" + region, binning="log")
    h_eff_seed_pt[region] = fill_hist_ratio(h_sim_pt_matchedSeed[region], h_sim_pt[region], "eff_Seed_pt_" + region, binning="log")
    h_eff_wrt_seed_pt[region] = fill_hist_ratio(h_sim_to_reco_match_pt[region], h_sim_pt_matchedSeed[region], "eff_wrt_seed_pt_" + region, binning = "log")

    h_eff_nrhits_byetaregion[region] = fill_hist_ratio(h_sim_nrhits_byetareg_matchedTrack[region], h_sim_nrhits_byetareg[region], "eff_nhits_" + region)

    h_eff_seed_nrhits_byetaregion[region] = fill_hist_ratio(h_sim_nrhits_byetareg_matchedSeed[region], h_sim_nrhits_byetareg[region], "eff_Seed_nhits_" + region)

    h_eff_nrhits_byetaregion_smallBrem[region] = fill_hist_ratio(h_sim_nrhits_byetareg_matchedTrack_smallBrem[region], h_sim_nrhits_byetareg_smallBrem[region], "eff_nhits_smallBrem" + region)

    h_eff_seed_nrhits_byetaregion_smallBrem[region] = fill_hist_ratio(h_sim_nrhits_byetareg_matchedSeed_smallBrem[region], h_sim_nrhits_byetareg_smallBrem[region], "eff_Seed_nhits_smallBrem" + region)
    


# -------------------- fake rates ---------------------------------
h_fakerate_eta = fill_hist_ratio(h_fake_eta, h_reco_eta,"fake_rate_eta")

h_fakerate_pt_barrel = fill_hist_ratio(h_fake_pt_barrel, h_reco_pt_barrel, "fake_rate_pt_barrel", binning = "log")
h_fakerate_pt_trans = fill_hist_ratio(h_fake_pt_trans, h_reco_pt_trans, "fake_rate_pt_trans", binning = "log")
h_fakerate_pt_endcap = fill_hist_ratio(h_fake_pt_endcap, h_reco_pt_endcap, "fake_rate_pt_endcap", binning = "log")

if args.is_gsf:
    outfile = "histograms/trackValHistograms" + args.ptmode + "GSF.root"
else:
    outfile = "histograms/trackValHistograms" + args.ptmode  + ".root"    


o = ROOT.TFile(outfile,"recreate")

h_reco_eta.Write()
h_fake_eta.Write()
h_sim_eta.Write()
h_sim_eta_trackerOnly.Write()
h_sim_eta_ecalOnly.Write()

h_sim_to_reco_match_eta.Write()
h_sim_eta_matchedSeed.Write()
h_sim_eta_matchedSeed_trackerOnly.Write()
h_sim_eta_matchedSeed_ecalOnly.Write

h_sim_nrhits.Write()
h_sim_nrhits_matchedTrack.Write()
h_sim_nrhits_matchedSeed.Write()
h_sim_nrhits_matchedSeed_trackerOnly.Write()
h_sim_nrhits_matchedSeed_ecalOnly.Write()
h_sim_nrhits_smallBrem.Write()
h_sim_nrhits_matchedTrack_smallBrem.Write()
h_sim_nrhits_matchedSeed_smallBrem.Write()
h_sim_nrhits_matchedSeed_smallBrem_trackerOnly.Write()
h_sim_nrhits_matchedSeed_smallBrem_ecalOnly.Write()

h_reco_pt_barrel.Write()
h_reco_pt_endcap.Write()
h_reco_pt_trans.Write()

h_fake_pt_barrel.Write()
h_fake_pt_endcap.Write()
h_fake_pt_trans.Write()

#------------efficiencies and fake rates wrt pt
for eta_region in eta_regions:
    h_eff_pt[eta_region].Write()
    h_eff_seed_pt[eta_region].Write()
    h_eff_wrt_seed_pt[eta_region].Write()

    h_sim_pt_matchedSeed[eta_region].Write()
    h_sim_to_reco_match_pt[eta_region].Write()
    h_sim_pt[eta_region].Write()

    h_eff_nrhits_byetaregion[eta_region].Write()
    h_eff_seed_nrhits_byetaregion[eta_region].Write()

    h_eff_nrhits_byetaregion_smallBrem[eta_region].Write()
    h_eff_seed_nrhits_byetaregion_smallBrem[eta_region].Write()

    h_sim_nrhits_byetareg[eta_region].Write()
    h_sim_nrhits_byetareg_matchedSeed[eta_region].Write()
    h_sim_nrhits_byetareg_matchedTrack[eta_region].Write()

h_fakerate_pt_barrel.Write()
h_fakerate_pt_endcap.Write()
h_fakerate_pt_trans.Write()
#-----efficiencies and fake rates wrt eta-----------

h_fakerate_eta.Write()
h_eff_eta.Write()
h_eff_seed_eta.Write()
h_eff_seed_eta_trackerOnly.Write()
h_eff_seed_eta_ecalOnly.Write()
h_eff_wrt_seed_eta.Write()

h_eff_nrhits.Write()
h_eff_nrhits_smallBrem.Write()
h_eff_seed_nrhits.Write()
h_eff_seed_nrhits_smallBrem.Write()

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
