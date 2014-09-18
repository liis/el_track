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
            eff = h_numerator.GetBinContent(i)/h_denominator.GetBinContent(i)
            h_ratio.SetBinContent(i, eff )
            h_ratio.SetBinError(i, ROOT.sqrt( eff*(1-eff)/h_denominator.GetBinContent(i) ) ) # uncertainty on efficiency

    return h_ratio

def fill_hist_ratio_poisson(h_numerator, h_denominator, outhistname, binning = "const"):
    """
    efficiency histogram with Poisson errorbars
    """
    #h_ratio_gr = ROOT.TGraphAsymmErrors(h_numerator.GetNbinsX())
    h_ratio_gr = ROOT.TGraphAsymmErrors(h_numerator)

    h_ratio_gr.BayesDivide(h_numerator, h_denominator)
    h_ratio_gr.SetName(outhistname)

    return h_ratio_gr


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

def initialize_histogram_2D(var, bin_reg_1, bin_reg_2 ):
    if len(bin_reg_1) == 3: # assume const. binning
        hist = ROOT.TH2F(var, var, bin_reg_1[0], bin_reg_1[1], bin_reg_1[2], bin_reg_2[0], bin_reg_2[1], bin_reg_2[2])
    elif len(bin_reg_1) == 2 and len(bin_reg_1[1]) > 2: # assume varying bin size
        hist = ROOT.TH2F(var, var, bin_reg_1[0], bin_reg_1[1], bin_reg_2[0], bin_reg_2[1])
    else:
        "Wrongly initialized histograms, please fix bin_reg and bin_type to allowed values"
        sys.exit()
    return hist


def initialize_histograms(vars, bin_reg, bin_reg_2 = "", hist_in_regions = [], dim = "1D" ):
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
                if dim == "2D":
                    histograms[var + "_" + region] = initialize_histogram_2D(var + "_" + region, bin_reg, bin_reg_2)
                else:
                    histograms[var + "_" + region] = initialize_histogram(var + "_" + region, bin_reg)


        else:
            if dim == "2D":
                histograms[var] = initialize_histogram_2D(var, bin_reg, vars[var])
            else:
                histograms[var] = initialize_histogram(var, bin_reg)
    return histograms
