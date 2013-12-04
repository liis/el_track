import ROOT

def fill_hist_ratio(h_numerator, h_denominator, const_bin_size=True):
    """
    bin-by-bin division of the input histograms
    """
    nbin = h_numerator.GetSize() #undrflow and overflow included
    xmin = h_numerator.GetXaxis().GetXmin()
    xmax = h_numerator.GetXaxis().GetXmax()
    print "xmax = ", xmax

    if const_bin_size:
        h_ratio = ROOT.TH1F("bla","bla",nbin-2,xmin,xmax)
    
    for i in range(1, nbin-1):
        if h_denominator.GetBinContent(i) == 0:
            h_ratio.SetBinContent(i,0)
        else:
            h_ratio.SetBinContent(i, h_numerator.GetBinContent(i)/h_denominator.GetBinContent(i) )

    return h_ratio

