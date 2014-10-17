import sys, os
import ROOT

import tdrstyle
tdrstyle.tdrstyle()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--infile', dest='infile', required=True) #directory for input histograms
args = parser.parse_args()

timing_maxCand = ROOT.TH1F("tn","tn", 7, 0, 7)
timing_maxCand.GetXaxis().SetTitle("max. cand")
timing_maxCand.GetYaxis().SetTitle("time (ns)")
timing_maxCand.GetYaxis().SetTitleOffset(1.6)
timing_maxCand.SetMinimum(0)

timing_nSigma = ROOT.TH1F("tns", "tns", 5, 0, 5)
timing_nSigma.GetXaxis().SetTitle("nSigma")
timing_nSigma.GetYaxis().SetTitle("time (ns)")
timing_nSigma.SetMinimum(0)

timing_chi2 = ROOT.TH1F("tnch", "tnch", 5, 0, 5)
timing_chi2.GetXaxis().SetTitle("max. chi2")
timing_chi2.GetYaxis().SetTitle("time (ns)")
timing_chi2.SetMinimum(0)

ichi2 = 0
chi2 = []

print "Opening timing file: " + args.infile
tf = ROOT.TFile(args.infile)
th = tf.Get("timing")

for ibin in range(1, th.GetNbinsX()+1 ):
    #    print "time at " + th.GetXaxis().GetBinLabel(ibin ) + " = " + str(th.GetBinContent(ibin))
    
    if  th.GetXaxis().GetBinLabel(ibin ): # if axis label is set
        if (th.GetXaxis().GetBinLabel(ibin).split("maxChi2_")[1]).split("_")[0] == str(2000) and (th.GetXaxis().GetBinLabel(ibin).split("nSigma_")[1]).split("_")[0] == str(3): #max cand
#            print "time at " + th.GetXaxis().GetBinLabel(ibin ) + " = " + str(th.GetBinContent(ibin))
            max_cand = (th.GetXaxis().GetBinLabel(ibin).split("maxCand_")[1]).split("_")[0]
            
            timing_maxCand.SetBinContent(int(max_cand), th.GetBinContent(ibin)*1000 )
            timing_maxCand.GetXaxis().SetBinLabel(int(max_cand), str(max_cand) )


        if (th.GetXaxis().GetBinLabel(ibin).split("maxChi2_")[1]).split("_")[0] == str(2000) and (th.GetXaxis().GetBinLabel(ibin).split("maxCand_")[1]).split("_")[0] == str(5): #nsigmna
#             print "time at " + th.GetXaxis().GetBinLabel(ibin ) + " = " + str(th.GetBinContent(ibin))

             nSigma = ((th.GetXaxis().GetBinLabel(ibin).split("nSigma_")[1]).split("_")[0])
    
             timing_nSigma.SetBinContent( int(nSigma), th.GetBinContent(ibin)*1000)
             timing_nSigma.GetXaxis().SetBinLabel( int(nSigma), str(nSigma))
             
    
        if (th.GetXaxis().GetBinLabel(ibin).split("nSigma_")[1]).split("_")[0] == str(3) and (th.GetXaxis().GetBinLabel(ibin).split("maxCand_")[1]).split("_")[0] == str(5): #nsigmna                                                 

            print "time at " + th.GetXaxis().GetBinLabel(ibin ) + " = " + str(th.GetBinContent(ibin))
            chi2_value = (th.GetXaxis().GetBinLabel(ibin).split("Chi2_")[1]).split("_")[0]

            chi2.append( [ int(chi2_value), th.GetBinContent(ibin)] ) # chi2 = [ [chi2_value, timing], [] ..], convert arg1 to int for sorting
            
chi2.sort()
print chi2

for ibin in range(0, len(chi2) ):
    timing_chi2.SetBinContent( ibin+1, chi2[ibin][1]*1000 )
    timing_chi2.GetXaxis().SetBinLabel( ibin+1, str(chi2[ibin][0]) )

c = ROOT.TCanvas("c1","c1", 800, 800)
timing_maxCand.Draw()            
c.SaveAs("timing/timing_maxCand.pdf")
c.Close()

c = ROOT.TCanvas("c1","c1", 800, 800)
timing_nSigma.Draw()
c.SaveAs("timing/timing_nSigma.pdf")
c.Close()

c = ROOT.TCanvas("c1","c1", 800, 800)
timing_chi2.Draw()
c.SaveAs("timing/timing_chi2.pdf")
c.Close()

