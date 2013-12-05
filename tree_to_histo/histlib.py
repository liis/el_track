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
        if h_denominator.GetBinContent(i) == 0:
            h_ratio.SetBinContent(i,0)
        else:
            h_ratio.SetBinContent(i, h_numerator.GetBinContent(i)/h_denominator.GetBinContent(i) )

    return h_ratio

