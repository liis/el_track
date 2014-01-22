import sys, os
import tdrstyle
tdrstyle.tdrstyle()
import ROOT
from plot_fun import infilenames_eta, infilenames_pt, infilenames_eta_gsf, infilenames_pt_gsf, load_input_files, draw_efficiency_histograms, draw_legend, draw_and_save_eff

indir = "$WORKING_DIR/tree_to_histo/histograms/" #location of input histograms -->

infile_eta = ROOT.TFile(indir + "trackValHistogramsPt100GSF.root")


h_nhit_tot = infile_eta.Get("sim_nrhits_matchedSeed_smallBrem")
h_nhit_ecalOnly = infile_eta.Get("sim_nrhits_matchedSeed_smallBrem_ecalOnly")
h_nhit_trackerOnly = infile_eta.Get("sim_nrhits_matchedSeed_smallBrem_trackerOnly")

hists = {"total": h_nhit_tot, "ECAL only": h_nhit_ecalOnly, "tracker only": h_nhit_trackerOnly}

h_nhit_tot.SetLineColor(ROOT.kBlack)
h_nhit_ecalOnly.SetLineColor(ROOT.kBlue)
h_nhit_trackerOnly.SetLineColor(ROOT.kRed)

c = ROOT.TCanvas("c","c", 600, 600)
c.SetGrid()
c.SetLogy()
xtitle = "nr. sim hits"
h_nhit_tot.GetXaxis().SetTitle(xtitle)

h_nhit_tot.Draw()
h_nhit_ecalOnly.Draw("same")
h_nhit_trackerOnly.Draw("same")

leg = draw_legend(hists, pos = "up_right")
leg.Draw()

c.SaveAs("$WORKING_DIR/plot/out_plots/comp_ecal_tracker_seed_smallBrem.pdf")
c.SaveAs("$WORKING_DIR/plot/out_plots/comp_ecal_tracker_seed_smallBrem.png")
#c.Close()
