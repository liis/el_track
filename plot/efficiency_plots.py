import sys
import os

import tdrstyle
tdrstyle.tdrstyle()
import ROOT

indir = "$WORKING_DIR/plot/input/elKF/" #location of input histograms

infiles_eff_pt = {"barrel": ["elFlatPtBarrelForEff.root"],
                  "endcap": ["elFlatPtEndcapForEff.root"], 
                  "trans": ["elFlatPtTransitionForEff.root"] 
                  }

infiles_fake_pt = {"barrel": ["elFlatPtBarrelForFake.root"],
                   "endcap": ["elFlatPtEndcapForFake.root"],
                   "trans": ["elFlatPtTransitionForFake.root"]
                   }

style = {"barrel": [1, 24], "endcap": [ROOT.kBlue, 25], "trans": [ROOT.kRed, 26] } #[color, markerStyle]


infiles_eff = {}
infiles_fake = {}
eff_hists = {}
fake_hists = {}
for region in infiles_eff_pt:
    inpath_eff = indir + infiles_eff_pt[region][0]
    inpath_fake = indir + infiles_fake_pt[region][0]
    print "loading efficiency histograms from: " + inpath_eff
    print "loading fake rate histograms from: " + inpath_fake 
    infiles_eff[region] = ROOT.TFile(inpath_eff) # need to write in dictionary, otherwise doesnt work
    infiles_fake[region] = ROOT.TFile(inpath_fake)
    
    eff_hists[region] = infiles_eff[region].Get("DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits/efficPt")
    fake_hists[region] = infiles_fake[region].Get("DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits/fakeratePt")


#------------ plot efficiencies----------------------
c_effpt = ROOT.TCanvas("c_eff","")
c_effpt.SetLogx()
c_effpt.SetGrid()

first = True
for region, hist in eff_hists.items():
    hist.SetLineColor( style[region][0])
    hist.SetMarkerStyle( style[region][1])
    hist.SetMarkerColor( style[region][0])
    if first:
        
        hist.GetXaxis().SetTitle("p_{T} (GeV)")
        hist.SetMaximum(1)
        hist.SetMinimum(0.)
        first = False
        hist.Draw()
    else:
        hist.Draw("same")

c_effpt.SaveAs("$WORKING_DIR/plot/out_plots/eff_el_pt.pdf")

#-----------------plot fake rates------------------------
c_fakept = ROOT.TCanvas("c_fake","")
c_fakept.SetLogx()
c_fakept.SetGrid()

first = True
for region, hist in fake_hists.items():
    hist.SetLineColor( style[region][0])
    hist.SetMarkerStyle( style[region][1])
    hist.SetMarkerColor( style[region][0])
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
                                                        
