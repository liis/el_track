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


def fill_hists_by_eta_regions(var_eta, var_to_hist, varname, hist_dict):
    """
    Fill histograms of var_to_hist for trans, endcap and barrel regions. pre-initialize three corresponding histograms in a dictionary
    """
    if abs(var_eta) < 0.9:
        hist_dict[varname + "_barrel"].Fill(var_to_hist)
    elif abs(var_eta) < 1.4:
        hist_dict[varname + "_trans"].Fill(var_to_hist)
    elif abs(var_eta) < 2.5:
        hist_dict[varname + "_endcap"].Fill(var_to_hist)
    
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

def initialize_histogram(var, bin_reg ):
    if len(bin_reg) == 3: # assume const. binning
        hist = ROOT.TH1F(var, var, bin_reg[0], bin_reg[1], bin_reg[2])
    elif len(bin_reg) == 2 and len(bin_reg[1]) > 2: # assume varying bin size
        hist = ROOT.TH1F(var, var, bin_reg[0], bin_reg[1])
    else:
        "Wrongly initialized histograms, please fix bin_reg and bin_type to allowed values"
        sys.exit()
    return hist


def initialize_histograms(vars, bin_reg, hist_in_regions = []):
    """
    Initialize histograms for a set of variables
    var -- list of histogram variables [var1, var2, ...]
    bin_reg -- (nbin, minbin, maxbin) for const or (nbin, xbinspt) for varying (e.g log)
    hist_in_regions -- define multiple histograms in different regions for the same variable [reg1, reg2, reg3]
    """

    histograms = {}
    for var in vars:
        if len(hist_in_regions) != 0:
            for region in hist_in_regions:
                histograms[var + "_" + region] = initialize_histogram(var + "_" + region, bin_reg)

        else:
            histograms[var] = initialize_histogram(var, bin_reg)
    return histograms



def initialize_hist_byEtaReg(histname, nbin, minbin, maxbin):
    hist_dict = {}
    for eta_region in eta_regions:
        hist_dict[eta_region] = ROOT.TH1I(histname + "_" + eta_region, histname + "_" + eta_region, nbin, minbin, maxbin)

    return hist_dict

