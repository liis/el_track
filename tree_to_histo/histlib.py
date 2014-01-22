import ROOT
import sys
from array import array

def log_binning(nbin,xmin,xmax):
    """
    convert histogram input to a list of logarithmically varying binsizes
    """
    logxmin = ROOT.TMath.log10(xmin)
    logxmax = ROOT.TMath.log10(xmax)
    binwidth = (logxmax - logxmin)/nbin
    xbins = []
    xbins.append(xmin)
    for i in range(1,nbin+1):
        xbins.append(xmin + ROOT.TMath.Power(10, logxmin + i*binwidth) )
    return xbins


def fill_hists_by_eta_regions(var_eta, var_to_hist, hist_dict):
    """
    Fill histograms of var_to_hist for trans, endcap and barrel regions. pre-initialize three corresponding histograms in a dictionary
    """
    if abs(var_eta) < 0.9:
        hist_dict["barrel"].Fill(var_to_hist)
    elif abs(var_eta) < 1.4:
        hist_dict["trans"].Fill(var_to_hist)
    elif abs(var_eta) < 2.5:
        hist_dict["endcap"].Fill(var_to_hist)
    
    return 0

def fill_hist_ratio(h_numerator, h_denominator, outhistname, binning = "const"):
    """
    bin-by-bin division of the input histograms
    Set const_bin_size = False for varing logarithmic bin size
    """
    nbin = h_numerator.GetSize() #undrflow and overflow included
    xmin = h_numerator.GetXaxis().GetXmin()
    xmax = h_numerator.GetXaxis().GetXmax()

    if binning == "const":
        h_ratio = ROOT.TH1F(outhistname,outhistname,nbin-2,xmin,xmax)
    elif binning == "log":
        h_ratio = ROOT.TH1F(outhistname, outhistname, nbin-2, array('d',log_binning(nbin-2,xmin,xmax)) )
    else:
        print "Error:fill_hist_ratio -- Implemented values for binning are: *const* or *log* "
        sys.exit()
        

    for i in range(1, nbin-1):
#        print "numerator = " + str(h_numerator.GetBinContent(i)) + ", denominator = " + str(h_denominator.GetBinContent(i))
        if h_denominator.GetBinContent(i) == 0:
            h_ratio.SetBinContent(i,0)
        else:
            h_ratio.SetBinContent(i, h_numerator.GetBinContent(i)/h_denominator.GetBinContent(i) )

#    h_ratio.Write()

    return h_ratio

eta_regions = ["barrel", "trans", "endcap"]

def initialize_hist_byEtaReg(histname, nbin, minbin, maxbin):
    hist_dict = {}
    for eta_region in eta_regions:
        hist_dict[eta_region] = ROOT.TH1I(histname + "_" + eta_region, histname + "_" + eta_region, nbin, minbin, maxbin)

    return hist_dict

