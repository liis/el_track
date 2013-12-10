import sys, os
import tdrstyle
tdrstyle.tdrstyle()
import ROOT
from plot_fun import infilenames_eta, infilenames_pt, load_input_files, get_hit_efficiency_hist, draw_efficiency_histograms

indir = "$WORKING_DIR/tree_to_histo/histograms/"

infiles_eta = load_input_files(indir, infilenames_eta)
infiles_pt = load_input_files(indir, infilenames_pt)

hist_hit_eff = {}
hist_hit_eff[ "eta_065"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 65, nrhits = 14)
hist_hit_eff["eta_075"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 75, nrhits = 14)
hist_hit_eff["eta_085"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 85, nrhits = 14)
hist_hit_eff["eta_095"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 95, nrhits = 14)

hist_hit_eff["pt_065"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 65, nrhits = 14)
hist_hit_eff["pt_075"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 75, nrhits = 14)
hist_hit_eff["pt_085"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 85, nrhits = 14)
hist_hit_eff["pt_095"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 95, nrhits = 14)

#---------------plots--------------------

for quality, hist_dictionary in hist_hit_eff.iteritems():
    for region, hists in hist_dictionary.iteritems():
        print "drawing histograms in region " + region

        c_eff = ROOT.TCanvas("c_eff_" + region + "_" +quality + "_" + region, "")
        c_eff.SetGrid()
        if region[:6] == "FlatPt":
            c_eff.SetLogx()

        draw_efficiency_histograms(hists,region)
        
        c_eff.SaveAs("$WORKING_DIR/plot/out_plots/hiteff_" + quality + "_" + region + ".pdf")
