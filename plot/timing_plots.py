import sys, os
import ROOT

import tdrstyle
tdrstyle.tdrstyle()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--infile', dest='infile', required=True) #directory for input histograms
args = parser.parse_args()

timing_maxCand = ROOT.TH1F("tn","tn", 8, 0, 8)

print "Opening timing file: " + args.infile
tf = ROOT.TFile(args.infile)
th = tf.Get("timing")

for ibin in range(1, th.GetNbinsX()+1 ):
    #    print "time at " + th.GetXaxis().GetBinLabel(ibin ) + " = " + str(th.GetBinContent(ibin))
    
    if  th.GetXaxis().GetBinLabel(ibin ): # if axis label is set
        if (th.GetXaxis().GetBinLabel(ibin).split("maxChi2_")[1]).split("_")[0] == str(2000) and (th.GetXaxis().GetBinLabel(ibin).split("nSigma_")[1]).split("_")[0] == str(3): #max cand
            print "time at " + th.GetXaxis().GetBinLabel(ibin ) + " = " + str(th.GetBinContent(ibin))

            max_cand = (th.GetXaxis().GetBinLabel(ibin).split("maxCand_")[1]).split("_")[0]
            
            timing_maxCand.SetBinContent(ibin, th.GetBinContent(ibin) )
            timing_maxCand.GetXaxis().SetBinLabel(ibin, str(max_cand) )

timing_maxCand.Draw()
