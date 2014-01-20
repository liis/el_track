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

fake_hists_eta = {}
for ptregion, infile_eta in infiles_eta.iteritems():
    histname_eff_eta = "eff_eta"
    histname_eff_seed_eta = "eff_seed_eta"

    histname_fake_eta = "fake_rate_eta"

    eff_hists_eta[ptregion] = infile_eta.Get(histname_eff_eta)
    eff_seed_hists_eta[ptregion] = infile_eta.Get(histname_eff_seed_eta)

    fake_hists_eta[ptregion] = infile_eta.Get(histname_fake_eta)

print "Getting histograms from files "
eff_hists_pt = {}
eff_seed_hists_pt = {}

fake_hists_pt = {}
print infiles_pt.values()
for infile_pt in infiles_pt.values():
    histname_eff_pt = "eff_pt"
    histname_eff_seed_pt = "eff_Seed_pt"

    histname_fake_pt = "fake_rate_pt"

    eff_hists_pt["barrel"] = infile_pt.Get(histname_eff_pt + "_barrel")
    eff_hists_pt["endcap"] = infile_pt.Get(histname_eff_pt + "_endcap")
    eff_hists_pt["trans"] = infile_pt.Get(histname_eff_pt + "_trans")

    eff_seed_hists_pt["barrel"] = infile_pt.Get(histname_eff_seed_pt + "_barrel")
    eff_seed_hists_pt["endcap"] = infile_pt.Get(histname_eff_seed_pt + "_endcap")
    eff_seed_hists_pt["trans"] = infile_pt.Get(histname_eff_seed_pt + "_trans")

    fake_hists_pt["barrel"] = infile_pt.Get(histname_fake_pt + "_barrel")
    fake_hists_pt["endcap"] = infile_pt.Get(histname_fake_pt + "_endcap")
    fake_hists_pt["trans"] = infile_pt.Get(histname_fake_pt + "_trans")

 

style = {"pt": {"barrel": [1, 24], "endcap": [ROOT.kRed, 25], "trans": [ROOT.kBlue, 26] },
         "eta": {"Pt1": [1, 24], "Pt10": [ROOT.kBlue, 25], "Pt100": [ROOT.kRed, 26] }}

print "Plotting efficiencies and fake rates"
draw_and_save_eff(eff_hists_eta, "eta", "eff", is_gsf = is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_hists_pt, "pt", "eff", is_gsf = is_gsf, leg_pos = "down_right")

draw_and_save_eff(eff_seed_hists_eta, "eta", "eff_seed", is_gsf = is_gsf, leg_pos="down_right")
draw_and_save_eff(eff_seed_hists_pt, "pt", "eff_seed", is_gsf = is_gsf, leg_pos="down_right")

draw_and_save_eff(fake_hists_eta, "eta", "fake", is_gsf = is_gsf, leg_pos = "up_right")
draw_and_save_eff(fake_hists_pt, "pt", "fake", is_gsf = is_gsf, leg_pos = "up_right")

