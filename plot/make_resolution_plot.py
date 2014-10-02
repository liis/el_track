import sys, os
import tdrstyle

tdrstyle.tdrstyle()
import ROOT
from plot_fun import infilenames_eta, infilenames_pt, infilenames_eta_gsf, infilenames_pt_gsf, load_input_files, draw_efficiency_histograms, draw_legend, draw_and_save_eff, draw_resolution

indir = "$WORKING_DIR/tree_to_histo/histograms/13TeV_v1/" #location of input histograms -->

#infile_eta = ROOT.TFile(indir + "trackValHistogramsPt100GSF.root")
infile_eta = ROOT.TFile(indir + "efficiencyHistograms_Pt10_maxCand_5_maxChi2_50_nSigma_3.root")

res_dxy_vs_eta = infile_eta.Get("resolutions/res_dxy_vs_eta")

#res_dxy_vs_eta.Draw("lego")

res = draw_resolution(res_dxy_vs_eta, "res_dxy_vs_eta")
res.Draw()

"""
xbins =  res_dxy_vs_eta.GetNbinsX() 
ybins = res_dxy_vs_eta.GetNbinsY()

res_dxy_vs_eta_out = ROOT.TH1F("test","test", xbins-1, -2.5, 2.5)

print "xbins = " + str(xbins)
print "ybins = " + str(ybins)



temp = ROOT.TH1F("res_test","res_test", ybins, res_dxy_vs_eta.GetBinLowEdge(1), res_dxy_vs_eta.GetBinLowEdge(xbins+1) )

print res_dxy_vs_eta.GetXaxis().GetXmin()
print res_dxy_vs_eta.GetXaxis().GetXmax()

#xBins = th2->GetNbinsX();
for i in range(1, xbins+1 ):
    temp = res_dxy_vs_eta.ProjectionY("proj",i, i+1)
    temp_gaus = ROOT.TF1("temp_gaus","gaus", temp.GetMean()-1*temp.GetRMS(), temp.GetMean()+1*temp.GetRMS() )#, -0.05, 0.05)

    r = temp.Fit(temp_gaus, "SRML Q")
    mean = r.Parameter(1)
    sigma = r.Parameter(2)

    res_dxy_vs_eta_out.SetBinContent(i, sigma)


res_dxy_vs_eta_out.Draw()
"""

#    if i == 1:
#        rms = temp.GetRMS()
#        temp.Draw()
#        print "rms = " + str(rms)
#        temp_gaus.Draw("same")
        
#        r = temp.Fit(temp_gaus, "SRML Q")

#        tmpmean = r.Parameter(1)
#        tmpsigma = r.Parameter(2)

#        print "sigma = " + str(tmpsigma)
#        temp_gaus.Draw()

#        xmin = mean - sigma*5
#        xmax = mean + sigma*5

#        meanRoo = ROOT.RooRealVar("mean", "mean of gaussian", tmpmean, xmin, xmax)

#        x = ROOT.RooRealVar("x","x", xmin, xmax)

        

#c = ROOT.TCanvas("c","c", 600, 600)   

"""
xvar = "pt"
extra = "_trans"

h_nhit_tot = infile_eta.Get("sim_" + xvar + "_matchedSeed" + extra)
h_nhit_ecalOnly = infile_eta.Get("sim_" + xvar + "_matchedSeed_ecalOnly" + extra)
h_nhit_trackerOnly = infile_eta.Get("sim_" + xvar + "_matchedSeed_trackerOnly" + extra)

hists = {"total": h_nhit_tot, "ECAL only": h_nhit_ecalOnly, "tracker only": h_nhit_trackerOnly}

h_nhit_tot.SetLineColor(ROOT.kBlack)
h_nhit_ecalOnly.SetLineColor(ROOT.kBlue)
h_nhit_trackerOnly.SetLineColor(ROOT.kRed)
"""


"""
c = ROOT.TCanvas("c","c", 600, 600)
c.SetGrid()
c.SetLogy()
c.SetLogx()
xtitle = xvar
h_nhit_tot.GetXaxis().SetTitle(xtitle)
h_nhit_tot.SetMinimum(1)

h_nhit_tot.Draw()
h_nhit_ecalOnly.Draw("same")
h_nhit_trackerOnly.Draw("same")

leg = draw_legend(hists, pos = "up_left")
leg.Draw()

c.SaveAs("$WORKING_DIR/plot/out_plots/comp_" + xvar + "_ecal_tracker_seed_smallBrem" + extra + ".pdf")
c.SaveAs("$WORKING_DIR/plot/out_plots/comp_" + xvar + "_ecal_tracker_seed_smallBrem" + extra + ".png")
#c.Close()
"""
