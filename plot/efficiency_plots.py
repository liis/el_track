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
eff_hists_eta = {}
eff_seed_hists_eta = {}
eff_wrt_seed_hists_eta = {} #nr matched to reco track/nr matched to reco seed

eff_hists_nrhits = {}
eff_seed_hists_nrhits = {}

eff_hists_nrhits_smallBrem = {}
eff_seed_hists_nrhits_smallBrem = {}

fake_hists_eta = {}
for ptregion, infile_eta in infiles_eta.iteritems():
    histname_eff_eta = "eff_eta"
    histname_eff_wrt_seed_eta = "eff_wrt_seed_eta"
    histname_eff_seed_eta = "eff_seed_eta"

    histname_eff_nrhits = "eff_nrhits"
    histname_eff_seed_nrhits = "eff_seed_nrhits"
    histname_eff_nrhits_smallBrem = "eff_nrhits_smallBrem"
    histname_eff_seed_nrhits_smallBrem = "eff_seed_nrhits_smallBrem"

    histname_fake_eta = "fake_rate_eta"

    eff_hists_eta[ptregion] = infile_eta.Get(histname_eff_eta)
    eff_wrt_seed_hists_eta[ptregion] = infile_eta.Get(histname_eff_wrt_seed_eta)
    eff_seed_hists_eta[ptregion] = infile_eta.Get(histname_eff_seed_eta)

    eff_hists_nrhits[ptregion] = infile_eta.Get(histname_eff_nrhits)
    eff_seed_hists_nrhits[ptregion] = infile_eta.Get(histname_eff_seed_nrhits)

    eff_hists_nrhits_smallBrem[ptregion] = infile_eta.Get(histname_eff_nrhits_smallBrem)
    eff_seed_hists_nrhits_smallBrem[ptregion] = infile_eta.Get(histname_eff_seed_nrhits_smallBrem)

    fake_hists_eta[ptregion] = infile_eta.Get(histname_fake_eta)

print "Getting histograms from files "
eff_hists_pt = {}
eff_wrt_seed_hists_pt = {}
eff_seed_hists_pt = {}

fake_hists_pt = {}
etaregions = ["barrel", "trans", "endcap"]
print infiles_pt.values()
for infile_pt in infiles_pt.values():
    histname_eff_pt = "eff_pt"
    histname_eff_wrt_seed_pt = "eff_wrt_seed_pt"
    histname_eff_seed_pt = "eff_Seed_pt"

    histname_fake_pt = "fake_rate_pt"

    for etaregion in etaregions:
        eff_hists_pt[etaregion] = infile_pt.Get(histname_eff_pt + "_" + etaregion)
        eff_seed_hists_pt[etaregion] = infile_pt.Get(histname_eff_seed_pt + "_" + etaregion)
        eff_wrt_seed_hists_pt[etaregion] = infile_pt.Get(histname_eff_wrt_seed_pt + "_" + etaregion)

        fake_hists_pt[etaregion] = infile_pt.Get(histname_fake_pt + "_" + etaregion)

style = {"pt": {"barrel": [1, 24], "endcap": [ROOT.kRed, 25], "trans": [ROOT.kBlue, 26] },
         "eta": {"Pt1": [1, 24], "Pt10": [ROOT.kBlue, 25], "Pt100": [ROOT.kRed, 26] }}

print "Plotting efficiencies and fake rates"
draw_and_save_eff(eff_hists_eta, "eta", "eff", is_gsf = is_gsf, leg_pos="down_right", title= "Total efficiency")
draw_and_save_eff(eff_seed_hists_eta, "eta", "eff_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Seeding efficiency")
draw_and_save_eff(eff_wrt_seed_hists_eta, "eta", "eff_wrt_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Reconstruction eff wrt. seeding")

draw_and_save_eff(eff_hists_pt, "pt", "eff", is_gsf = is_gsf, leg_pos = "down_right", title="Total efficiency")
draw_and_save_eff(eff_seed_hists_pt, "pt", "eff_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Seeding efficiency")
draw_and_save_eff(eff_wrt_seed_hists_pt, "pt", "eff_wrt_seed", is_gsf = is_gsf, leg_pos="down_right", title = "Reconstruction eff wrt. seeding")

draw_and_save_eff(eff_hists_nrhits, "nrhits", "eff", is_gsf = is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_seed_hists_nrhits, "nrhits", "eff_seed", is_gsf = is_gsf, leg_pos="down_right")

draw_and_save_eff(eff_hists_nrhits_smallBrem, "nrhits", "eff_smallBrem", is_gsf = is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_seed_hists_nrhits_smallBrem, "nrhits", "eff_seed_smallBrem", is_gsf = is_gsf, leg_pos="down_right")

draw_and_save_eff(fake_hists_eta, "eta", "fake", is_gsf = is_gsf, leg_pos = "up_right")
draw_and_save_eff(fake_hists_pt, "pt", "fake", is_gsf = is_gsf, leg_pos = "up_right")

