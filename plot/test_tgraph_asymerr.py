import ROOT

indir = "$WORKING_DIR/tree_to_histo/histograms/" #location of input histograms
infile = ROOT.TFile(indir + "trackValHistogramsPt100GSF.root")

h_num = infile.Get("sim_nrhits_matchedTrack_smallBrem")
h_den = infile.Get("sim_nrhits_smallBrem")

gr = ROOT.TGraphAsymmErrors( h_num.GetNbinsX() )
#gr.SetName( str(i_draw ) )
gr.BayesDivide( h_num, h_den)
gr.GetXaxis().SetTitle("#eta")
gr.GetYaxis().SetTitle("Efficiency")
gr.SetTitle("")
gr.SetMarkerColor(ROOT.kBlue)
gr.SetLineColor(ROOT.kBlue)
gr.SetMarkerStyle(20)
gr.SetMarkerSize(1)

c = ROOT.TCanvas("c1", "c1", 600, 600)
c.SetGrid()

gr.Draw("AP")
