import sys, os
import tdrstyle
tdrstyle.tdrstyle()
import ROOT
from plot_fun import infilenames_eta, infilenames_pt, load_input_files, get_hit_efficiency_hist

indir = "$WORKING_DIR/tree_to_histo/histograms/"

infiles_eta = load_input_files(indir, infilenames_eta)
infiles_pt = load_input_files(indir, infilenames_pt)


hist_hit_eff_eta = get_hit_efficiency_hist(infiles_eta, "eta", quality = 95, nrhits = 13)


#---------------plots--------------------



for region, hists in hist_hit_eff_eta.iteritems():
    print "drawing histograms in region " + region
    n = 0
    c_eff_eta = ROOT.TCanvas("c_eff_eta_" + region, "")
    c_eff_eta.SetGrid()
    for hist in hists:
        print "drawing histogram at hit " + str(n)
        if n == 0:
            hist.SetMaximum(1.)
            hist.SetMinimum(0.)
            hist.Draw()
        else:
            hist.Draw("same")
        n = n + 1
    c_eff_eta.SaveAs("$WORKING_DIR/plot/out_plots/hiteff_eta_" + region + ".pdf")
