import sys, os
import tdrstyle
tdrstyle.tdrstyle()
import ROOT
from plot_fun import infilenames_eta, infilenames_pt, load_input_files

indir = "$WORKING_DIR/tree_to_histo/histograms/" #location of input histograms -->clean up!

infiles_eta = load_input_files(indir, infilenames_eta)
infiles_pt = load_input_files(indir, infilenames_pt)


fake_hists_pt = {}
for region in infiles_pt:
    histname_pt = "fake_rate_pt_" + str(region)
    print histname_pt
    fake_hists_pt[region] = infile_pt.Get(histname_pt)

eff_hists_pt = {}
for region in infiles_eff_pt:
    histname_pt = "eff_pt_" +str(region)
    eff_hists_pt[region] = infile_pt.Get(histname_pt)


eff_hists_eta = {}
for region in infiles_eta:
    inpath_eta = indir_eta + infiles_eta[region][0]
    eff_hists_eta[region] = infiles_eta[region].Get("eff_eta")
    fake_hists_eta[region] = infiles_eta[region].Get("fake_rate_eta")
"""
"""
#------------ plot efficiencies----------------------
style = {"pt": {"barrel": [1, 24], "endcap": [ROOT.kRed, 25], "trans": [ROOT.kBlue, 26] },
         "eta": {"Pt1": [1, 24], "Pt10": [ROOT.kBlue, 25], "Pt100": [ROOT.kRed, 26] }}

c_effpt = ROOT.TCanvas("c_eff","")
c_effpt.SetLogx()
c_effpt.SetGrid()

first = True
for region, hist in eff_hists_pt.items():
    hist.SetLineColor( style["pt"][region][0])
    hist.SetMarkerStyle( style["pt"][region][1])
    hist.SetMarkerColor( style["pt"][region][0])
    if first:       
        hist.GetXaxis().SetTitle("p_{T} (GeV)")
        hist.SetMaximum(1)
        hist.SetMinimum(0.)
        first = False
        hist.Draw()
    else:
        hist.Draw("same")

c_effpt.SaveAs("$WORKING_DIR/plot/out_plots/eff_el_pt.pdf")
c_effpt.Close()

c_effeta = ROOT.TCanvas("c_effeta","c_effeta")
c_effeta.SetGrid()

first = True
for region, hist in eff_hists_eta.items():
    hist.SetLineColor( style["eta"][region][0])
    hist.SetMarkerStyle( style["eta"][region][1])
    hist.SetMarkerColor( style["eta"][region][0])
                   
    if first:
        hist.GetXaxis().SetTitle("#eta")
        hist.SetMaximum(1)
        hist.SetMinimum(0.5)
        hist.Draw("")
        first = False
    else:
        hist.Draw("same")
    
c_effeta.SaveAs("$WORKING_DIR/plot/out_plots/eff_el_eta.pdf")
c_effeta.Close()

#-----------------plot fake rates------------------------
c_fakept = ROOT.TCanvas("c_fake","")
c_fakept.SetLogx()
c_fakept.SetGrid()

first = True
for region, hist in fake_hists_pt.items():
    hist.SetLineColor( style["pt"][region][0])
    hist.SetMarkerStyle( style["pt"][region][1])
    hist.SetMarkerColor( style["pt"][region][0])
    if first:
        hist.GetXaxis().SetTitle("p_{T} (GeV)")
#        hist.SetAxisRange(0.1, 150, "x") 
        hist.SetMaximum(0.3)
        hist.SetMinimum(0.)
        first = False
        hist.Draw()
    else:
        hist.Draw("same")
        
c_fakept.SaveAs("$WORKING_DIR/plot/out_plots/fake_el_pt.pdf")
                                                        
c_fake_eta = ROOT.TCanvas("c_fake_eta","")
c_fake_eta.SetGrid()

first = True
for region, hist in fake_hists_eta.items():
    hist.SetLineColor( style["eta"][region][0])
    hist.SetMarkerColor( style["eta"][region][0])
    hist.SetMarkerStyle( style["eta"][region][1])
    if first:
        hist.GetXaxis().SetTitle("#eta")
        hist.SetMaximum(0.3)
        hist.SetMinimum(0.)
        hist.Draw()
        first = False
    else:
        hist.Draw("same")
c_fake_eta.SaveAs("$WORKING_DIR/plot/out_plots/fake_el_eta.pdf")
"""
