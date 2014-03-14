import sys, os
import tdrstyle
tdrstyle.tdrstyle()
import ROOT
from plot_fun import infilenames_eta, infilenames_pt, infilenames_eta_gsf, infilenames_pt_gsf, load_input_files, draw_efficiency_histograms, draw_legend, draw_and_save_eff

indir = "$WORKING_DIR/tree_to_histo/histograms/" #location of input histograms -->clean up!

is_gsf = True

if is_gsf:
    infiles_eta = load_input_files(indir, infilenames_eta_gsf)
    infiles_pt = load_input_files(indir, infilenames_pt_gsf)
else:
    infiles_eta = load_input_files(indir, infilenames_eta)
    infiles_pt = load_input_files(indir, infilenames_pt)


print "Opening input files from " + indir

##-------------- Group histograms wrt. Pt regions (Pt10, Pt100)--------------
eff_eta = {}
eff_seed_eta = {}
eff_seed_charge_eta = {}
eff_wrt_seed_eta = {} #nr matched to reco track/nr matched to reco seed

eff_eta_smallBrem = {}
eff_seed_eta_smallBrem = {}
eff_wrt_seed_eta_smallBrem = {}

eff_seed_eta_tracker = {}
eff_seed_eta_ecal = {}

eff_nrhits = {}
eff_wrt_seed_nrhits = {}
eff_seed_nrhits = {}

eff_nrhits_smallBrem = {}
eff_seed_nrhits_smallBrem = {}
eff_wrt_seed_nrhits_smallBrem = {}

fake_eta = {}

pt_pull = {}
pt_res = {}
for ptregion, infile_eta in infiles_eta.iteritems():
    histname_eff_eta = "eff_eta"
    histname_eff_wrt_seed_eta = "eff_wrt_seed_eta"
    histname_eff_seed_eta = "eff_seed_eta"
    histname_eff_seed_charge_eta = "eff_seed_charge_eta"

    histname_eff_nrhits = "eff_nrhits"
    histname_eff_seed_nrhits = "eff_seed_nrhits"
    histname_eff_nrhits_smallBrem = "eff_nrhits_smallBrem"
    histname_eff_seed_nrhits_smallBrem = "eff_seed_nrhits_smallBrem"

    histname_fake_eta = "fake_rate_eta"
    histname_pt_pull = "pt_pull"
    histname_pt_res = "pt_res"

    eff_eta[ptregion] = infile_eta.Get(histname_eff_eta)
    eff_wrt_seed_eta[ptregion] = infile_eta.Get(histname_eff_wrt_seed_eta)
    eff_seed_eta[ptregion] = infile_eta.Get(histname_eff_seed_eta)
    eff_seed_charge_eta[ptregion] = infile_eta.Get(histname_eff_seed_charge_eta)
    eff_eta_smallBrem[ptregion] = infile_eta.Get("eff_eta_smallBrem")
    eff_wrt_seed_eta_smallBrem[ptregion] = infile_eta.Get("eff_wrt_seed_eta_smallBrem")
    eff_seed_eta_smallBrem[ptregion] = infile_eta.Get("eff_seed_eta_smallBrem")

    eff_seed_eta_tracker[ptregion] = infile_eta.Get("eff_seed_eta_tracker")
    eff_seed_eta_ecal[ptregion] = infile_eta.Get("eff_seed_eta_ecal")

    eff_nrhits[ptregion] = infile_eta.Get(histname_eff_nrhits)
    eff_seed_nrhits[ptregion] = infile_eta.Get(histname_eff_seed_nrhits)
    eff_wrt_seed_nrhits[ptregion] = infile_eta.Get("eff_wrt_seed_nrhits")

    eff_nrhits_smallBrem[ptregion] = infile_eta.Get(histname_eff_nrhits_smallBrem)
    eff_seed_nrhits_smallBrem[ptregion] = infile_eta.Get(histname_eff_seed_nrhits_smallBrem)
    eff_wrt_seed_nrhits_smallBrem[ptregion] = infile_eta.Get("eff_wrt_seed_nrhits_smallBrem")

    fake_eta[ptregion] = infile_eta.Get(histname_fake_eta)
    pt_pull[ptregion] = infile_eta.Get(histname_pt_pull)
    pt_res[ptregion] = infile_eta.Get(histname_pt_res)

#----------Group histograms wrt. eta regions (barrel, trans, endcap)--------------
eff_pt = {}
eff_wrt_seed_pt = {}
eff_seed_pt = {}
eff_seed_charge_pt = {}

eff_pt_smallBrem = {}
eff_wrt_seed_pt_smallBrem = {}
eff_seed_pt_smallBrem = {}

eff_seed_pt_ecal = {}
eff_wrt_seed_pt_ecal = {}
eff_seed_pt_eta = {}

eff_seed_pt_tracker = {}
eff_wrt_seed_pt_tracker = {}

eff_nrhits_byeta = {}
eff_seed_nrhits_byeta = {}
eff_wrt_seed_nrhits_byeta = {}

eff_nrhits_smallBrem_byeta = {}
eff_seed_nrhits_smallBrem_byeta = {}
eff_wrt_seed_nrhits_smallBrem_byeta = {}

fake_pt = {}
etaregions = ["barrel", "trans", "endcap"]
#print infiles_pt.values()
for infile_pt in infiles_pt.values():
    histname_eff_pt = "eff_pt"
    histname_eff_wrt_seed_pt = "eff_wrt_seed_pt"
    histname_eff_seed_pt = "eff_seed_pt"
    histname_eff_seed_charge_pt = "eff_seed_charge_pt"

    histname_fake_pt = "fake_rate_pt"

    for etaregion in etaregions:
        eff_pt[etaregion] = infile_pt.Get(histname_eff_pt + "_" + etaregion)
        eff_seed_pt[etaregion] = infile_pt.Get(histname_eff_seed_pt + "_" + etaregion)
        eff_seed_charge_pt[etaregion] = infile_pt.Get(histname_eff_seed_charge_pt + "_" + etaregion)
        eff_wrt_seed_pt[etaregion] = infile_pt.Get(histname_eff_wrt_seed_pt + "_" + etaregion)

        eff_pt_smallBrem[etaregion] = infile_pt.Get("eff_pt_smallBrem_" + etaregion)
        eff_seed_pt_smallBrem[etaregion] = infile_pt.Get("eff_seed_pt_smallBrem_" + etaregion)
        eff_wrt_seed_pt_smallBrem[etaregion] = infile_pt.Get("eff_wrt_seed_pt_smallBrem_" + etaregion)

        eff_seed_pt_tracker[etaregion] = infile_pt.Get("eff_seed_pt_tracker_" + etaregion)
        eff_wrt_seed_pt_tracker[etaregion] = infile_pt.Get("eff_wrt_seed_pt_tracker_" + etaregion)

        eff_seed_pt_ecal[etaregion]= infile_pt.Get("eff_seed_pt_ecal_" + etaregion)
        eff_wrt_seed_pt_ecal[etaregion] = infile_pt.Get("eff_wrt_seed_pt_ecal_" + etaregion)

        eff_nrhits_byeta[etaregion] = infile_pt.Get("eff_nhits_" + etaregion)
        eff_seed_nrhits_byeta[etaregion] = infile_pt.Get("eff_seed_nhits_" + etaregion)
        eff_wrt_seed_nrhits_byeta[etaregion] = infile_pt.Get("eff_wrt_seed_nhits_" + etaregion)

        eff_nrhits_smallBrem_byeta[etaregion] = infile_pt.Get("eff_nhits_smallBrem_" + etaregion)
        eff_seed_nrhits_smallBrem_byeta[etaregion] = infile_pt.Get("eff_seed_nhits_smallBrem_" + etaregion)
        eff_wrt_seed_nrhits_smallBrem_byeta[etaregion] = infile_pt.Get("eff_wrt_seed_nhits_smallBrem_" + etaregion)

        fake_pt[etaregion] = infile_pt.Get(histname_fake_pt + "_" + etaregion)

#----------------add custom groups if needed --------------------------------------


#-----------------------------------------------------------------------------------
style = {"pt": {"barrel": [1, 24], "endcap": [ROOT.kRed, 25], "trans": [ROOT.kBlue, 26] },
         "eta": {"Pt1": [1, 24], "Pt10": [ROOT.kBlue, 25], "Pt100": [ROOT.kRed, 26] }}

print "Plotting efficiencies and fake rates"
#---------------pt----------------------
draw_and_save_eff(eff_pt, "pt", "eff", is_gsf = is_gsf, leg_pos = "down_right", title="Total efficiency")
draw_and_save_eff(eff_seed_pt, "pt", "eff_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Seeding efficiency")
draw_and_save_eff(eff_wrt_seed_pt, "pt", "eff_wrt_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Reconstruction eff wrt. seeding")
draw_and_save_eff(eff_seed_charge_pt, "pt", "eff_seed_charge", is_gsf = is_gsf, leg_pos="down_right", title = "Seed charge efficiency")

draw_and_save_eff(eff_pt_smallBrem, "pt", "eff", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")
draw_and_save_eff(eff_seed_pt_smallBrem, "pt", "eff_seed",is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")
draw_and_save_eff(eff_wrt_seed_pt_smallBrem, "pt", "eff_wrt_seed", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")

draw_and_save_eff(eff_seed_pt_tracker,"pt", "eff_seed", is_gsf = is_gsf, leg_pos="down_right", label="trackerSeed")
draw_and_save_eff(eff_seed_pt_ecal, "pt", "eff_seed",is_gsf = is_gsf, leg_pos="down_right", label="ecalSeed")

#------------ nrhits -------------------
draw_and_save_eff(eff_nrhits, "nrhits", "eff", is_gsf=is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_seed_nrhits, "nrhits", "eff_seed", is_gsf = is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_wrt_seed_nrhits, "nrhits", "eff_wrt_seed", is_gsf = is_gsf, leg_pos = "down_right")

draw_and_save_eff(eff_nrhits_smallBrem, "nrhits", "eff", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")
draw_and_save_eff(eff_wrt_seed_nrhits_smallBrem, "nrhits", "eff_wrt_seed", is_gsf, label = "smallBrem", leg_pos = "down_right")
draw_and_save_eff(eff_seed_nrhits_smallBrem, "nrhits", "eff_seed", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")

draw_and_save_eff(eff_nrhits_byeta, "nrhits_byeta", "eff", is_gsf = is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_seed_nrhits_byeta, "nrhits_byeta", "eff_seed", is_gsf = is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_wrt_seed_nrhits_byeta, "nrhits_byeta", "eff_wrt_seed", is_gsf = is_gsf, leg_pos = "down_right")

draw_and_save_eff(eff_nrhits_smallBrem_byeta, "nrhits_byeta", "eff", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")
draw_and_save_eff(eff_wrt_seed_nrhits_smallBrem_byeta, "nrhits_byeta", "eff_wrt_seed", is_gsf, label = "smallBrem", leg_pos = "down_right")
draw_and_save_eff(eff_seed_nrhits_smallBrem_byeta, "nrhits_byeta", "eff_seed", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")


#----------- eta ---------------
draw_and_save_eff(eff_eta_smallBrem, "eta", "eff", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")
draw_and_save_eff(eff_seed_eta_smallBrem, "eta", "eff_seed", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")
draw_and_save_eff(eff_wrt_seed_eta_smallBrem, "eta", "eff_wrt_seed", is_gsf = is_gsf, label = "smallBrem", leg_pos="down_right")


draw_and_save_eff(eff_eta, "eta", "eff", is_gsf = is_gsf, leg_pos="down_right", title= "Total efficiency")
draw_and_save_eff(eff_seed_eta, "eta", "eff_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Seeding efficiency")
draw_and_save_eff(eff_wrt_seed_eta, "eta", "eff_wrt_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Reconstruction eff wrt. seeding")
draw_and_save_eff(eff_seed_charge_eta, "eta", "eff_seed_charge", is_gsf = is_gsf, leg_pos="down_right", title = "Seed charge identification efficiency")


draw_and_save_eff(eff_seed_eta_tracker, "eta", "eff_seed", is_gsf = is_gsf, leg_pos="down_right", label="trackerSeed")
draw_and_save_eff(eff_seed_eta_ecal, "eta", "eff_seed", is_gsf = is_gsf, leg_pos="down_right", label="ecalSeed")

#----------------- fake --------------
draw_and_save_eff(fake_eta, "eta", "fake", is_gsf = is_gsf, leg_pos = "up_right")
draw_and_save_eff(pt_pull, "pt_pull", "pull", is_gsf = is_gsf, leg_pos = "up_right")
draw_and_save_eff(pt_res, "pt_res", "res", is_gsf = is_gsf, leg_pos = "up_right")

draw_and_save_eff(fake_pt, "pt", "fake", is_gsf = is_gsf, leg_pos = "up_right")

