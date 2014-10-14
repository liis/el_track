import ROOT
from ROOT import RooFit
import math
from collections import OrderedDict as dict

infilenames_eta = {
    "Pt100": "trackValHistograms_trackValTree_Pt100_maxCand_5_maxChi2_2000_nSigma_3.root",
    "Pt10":  "trackValHistograms_trackValTree_Pt10_maxCand_5_maxChi2_2000_nSigma_3.root",
#    "Pt100": "trackValHistogramsPt100.root"
}

infilenames_eta_gsf = {
    "Pt100":   "trackValHistograms_trackValTree_Pt100_maxCand_5_maxChi2_2000_nSigma_3.root",
    "Pt10" :   "trackValHistograms_trackValTree_Pt10_maxCand_5_maxChi2_2000_nSigma_3.root",
#    "Pt10": "trackValHistogramsPt10GSF.root",
#    "Pt100": "trackValHistogramsPt100GSF.root"
}

infilenames_pt = {
#    "FlatPt": "trackValHistogramsFlatPt.root"
    "FlatPt": "trackValHistograms_trackValTree_FlatPt_maxCand_5_maxChi2_100_nSigma_3.root"
}

infilenames_pt_gsf = {
    "FlatPt": "trackValHistograms_trackValTree_FlatPt_maxCand_5_maxChi2_100_nSigma_3.root"
#    "FlatPt": "trackValHistogramsFlatPtGSF.root"
}

def load_input_files(indir, infilenames):
    """
    load input files. give dir path and dictionary of filenames as input
    """
    infiles = dict()

    for region, filename in infilenames.iteritems():
        inpath = indir + filename
        infiles[region] = ROOT.TFile(inpath)
        
    return infiles

def get_hit_efficiency_hist(infiles, var, quality, nrhits = 10):
    """
    read a dictionary of input root files in different regions ("Pt1", "Pt10", "Pt100", "FlatPt") and return a dictionary of histograms for each hit
    var -- eta or pt
    quality -- 65, 75, 85, 95 - corresponding to the pt efficiency required at each hit
    """

    hit_dir = "VariablesBySimhit/"
    hit_eff_hists = {}

    for region, histfile in infiles.iteritems():
        if region=="FlatPt":
            hit_eff_hists[region + "_barrel"] = []
            hit_eff_hists[region + "_trans"] = []
            hit_eff_hists[region + "_endcap"] = []
        else:
            hit_eff_hists[region] = [] # at each region save a vector of efficiencies for each hit

        for ihit in range(0, nrhits):
            if region=="FlatPt": #separate different eta regions
                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_barrel_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region+"_barrel"].append(histfile.Get(hist_name))

                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_trans_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region+"_trans"].append(histfile.Get(hist_name))

                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_endcap_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region+"_endcap"].append(histfile.Get(hist_name))
            else:
                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region].append( histfile.Get(hist_name) )



    return hit_eff_hists


colors = [ROOT.kBlack, ROOT.kRed, ROOT.kYellow, ROOT.kYellow-3, ROOT.kGreen, ROOT.kGreen+3, ROOT.kCyan, ROOT.kCyan+1, ROOT.kCyan+2, ROOT.kCyan+3, ROOT.kCyan+4, ROOT.kBlue, ROOT.kViolet-3, ROOT.kViolet+3]

def draw_efficiency_histograms(hists, xtitle = "none", ytitle = "none", ymax =  "", region="none", style="", logy=False, is_secondary=False):
    """
    hists -- dictionary of histograms
    plot a list of efficiency histogras
    secondary = True -- add additional histograms to existing ones
    """
    
    n = 0
    for hist in hists:
        hist.SetLineColor(colors[n])
        hist.SetMarkerColor(colors[n])
        hist.SetMarkerStyle(20);
        hist.SetMarkerSize(1);
        if n==0:
            if ymax == "":
                ymax=hist.GetMaximum()*2

            hist.SetMaximum(ymax)
            hist.SetMinimum(0.)
            if logy:
                hist.SetMinimum(hist.GetMinimum()*0.2)
            hist.GetXaxis().SetTitleOffset(1.)
            hist.GetYaxis().SetTitleOffset(1.2)
#            hist.SetTitle("blabla")
            hist.SetStats(False)

            if( region[:3]=="Pt1"):
                hist.GetXaxis().SetTitle("#eta")
            elif( xtitle=="p_{T}"):
                hist.SetAxisRange(1., 200, 'x')

            if not xtitle == "none":
                hist.GetXaxis().SetTitle(xtitle)
            if not ytitle == "none":
                hist.GetYaxis().SetTitle(ytitle)

#            print "style = " + style
            if style == "noerr":
                if is_secondary:
                    hist.Draw("histsame")
                else:
                    hist.Draw("hist")
            elif style == "":
                if is_secondary:
                    hist.Draw("epsame")
                else:
                    hist.Draw("ep")
            else:
                hist.SetMarkerStyle(style)
                if is_secondary:
                    hist.Draw("psame")
                else:
                    hist.Draw("p")
        else:
            if style == "noerr":
                hist.Draw("histsame")
            elif style == "":
                hist.Draw("epsame")
            else:
                hist.SetMarkerStyle(style)
                hist.Draw("psame")

        n = n + 1

def draw_legend(hists, pos = "down_right"):
    """
    hist - dictionary of process names and histograms
    """
    if pos == "down_right":
        leg = ROOT.TLegend(0.55,0.5,0.95,0.29);

    if pos == "up_right":
        leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.89);

    if pos == "up_left":
        leg = ROOT.TLegend(0.2, 0.7, 0.45, 0.89);

    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    for process,hist in hists.iteritems():
        leg.AddEntry(hist, process)

    return leg

def add_text_box(text=''):
    """ Add a CMS blurb to a plot """
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.SetTextAlign(11)
    return latex.DrawLatex(0.17, 0.96, text)


#def draw_and_save_res(hists, var, res_var, is_gsf, label = "", leg_pos = "up_right", title = "", ymax_res = 0.5, style=""):
#    c = ROOT.TCanvas("c","c", 600, 600)
#    c.SetGrid()

#    if var == "eta":
#        xtitle = "eta"

#    ytitle = "resolution in " + pretty_var[res+var]

res_map = {
    "dxy": "d_{0}",
    "dz": "z_{0}",
    "pt": "(res p_{T})/p_{T} [%]",
    "cotth": "cot(#theta)",
    "phi": "#phi",
    }
    
def draw_and_save_res(hists68, hists95, var, resvar, is_gsf, outdir, leg_pos = "up_right"):
    c = ROOT.TCanvas("c","c", 800, 800)
    c.SetGrid()
    if var == "eta":
        xtitle = "#eta"
    if var == "pt":
        xtitle = "p_{T}"

    ytitle = "Resolution in " + res_map[resvar]

    draw_efficiency_histograms(hists95.values(), xtitle, ytitle, style=4) #do 95 first, since the hist maximum is set here
    draw_efficiency_histograms(hists68.values(), xtitle, ytitle, style=20, is_secondary=True)

    leg_68 = draw_legend(hists68, pos = "up_right")
#    leg_68.SetHeader(legend_header)
    leg_68.Draw()

    GSFstr = ""
    if(is_gsf):
        GSFstr = "_GSF"

    c.SaveAs(outdir +"/resolution_" + var + "_" + resvar + GSFstr + ".png")
    c.Close()

        
def draw_and_save_eff(hists, var, eff_fake, is_gsf, label = "", leg_pos = "up_right", title = "", ymax_res=0.5, style=""):
    """
    hists - dictionary of process names and histograms
    var - xaxis variable
    eff_fake - "eff", "eff_seed", "eff_wrt_seed", "fake"
    label - "extra label"
    """

    c = ROOT.TCanvas("c","c", 600, 600)
    c.SetGrid()

    logy=False
    xtitle = ""
    if var == "pt":
        c.SetLogx()
        xtitle = "p_{T}"
    if var == "eta":
        xtitle = "#eta"
    if var[:6] == "nrhits":
        xtitle = "Number of sim. hits"

    ytitle = ""
    ymax = 1
    if eff_fake == "eff":
        ytitle = "Efficiency "
    if eff_fake == "eff_seed":
        ytitle = "Seeding efficiency"
    if eff_fake == "eff_wrt_seed":
        ytitle = "Reco wrt seeding efficiency"
    if eff_fake[:4] == "fake":
        ytitle = "Fake rate"
        if var=="pt":
            ymax = 0.25
        elif var=="eta":
            ymax = 0.25
        else:
            ymax = 1
#        if var == "pt":
#            ymax = ymax_res

    if eff_fake[:4] == "pull":
        ytitle = "pull"
    if eff_fake[:4] == "res":
        ytitle = label
        ymax = ymax_res
        logy=False # set nonzero y value
#        c.SetLogy

    if len(title)>0:
        ytitle=ytitle + " (" + title + ")"

#    if len(label) > 0:
#        ytitle = ytitle + " (" + label + ")"

    draw_efficiency_histograms(hists.values(), xtitle, ytitle, ymax, style=style, logy=logy)
    leg = draw_legend(hists, pos = leg_pos)
    leg.Draw()

    GSFstr = ""
    if(is_gsf):
        GSFstr = "_GSF"
    if logy:
        c.SetLogy()

#    c.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/" + eff_fake + "_" + var + "_" + label + GSFstr + ".pdf")
#    c.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/13TeV_011014_a/" + eff_fake + "_" + var + "_" + label + GSFstr + ".png")
#    c.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/13TeV_011014_eff/" + eff_fake + "_" + var + "_" + label + GSFstr + ".png")
    c.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/13TeV_011014_res/" + eff_fake + "_" + var + "_" + label + GSFstr + ".png")
    
    c.Close()


def draw_resolution(res_hist_2d, res_hist_name, do_control_fit_plots=False, mode="68"):
    """
    res_hist_2d -- histogram of eta/pt vs resolution variable
    res_hist_name -- name to be given to the output resolution histogram
    """
    print "draw control fits = " + str(do_control_fit_plots)

    nbinsx = res_hist_2d.GetNbinsX()
    nbinsy = res_hist_2d.GetNbinsY()

    temp_res = ROOT.TH1F("temp_res", "temp_res", nbinsy, res_hist_2d.GetBinLowEdge(1), res_hist_2d.GetBinLowEdge(nbinsx+1) )

    res_1d = ROOT.TH1F(res_hist_name,res_hist_name, nbinsx-1, res_hist_2d.GetXaxis().GetXmin(), res_hist_2d.GetXaxis().GetXmax() )

    for i in range(1, nbinsx+1):
        
        temp_res = res_hist_2d.ProjectionY("proj",i, i+1) #get resolution distribution at each bin
        
        isPtPlot = False
        #-- initial guess of fit boarders ---
        if "res_dxy" in res_hist_name:
          firstRangeLeft = -0.03
          firstRangeRight = 0.03
          tailSign = 0
        elif "res_dz" in res_hist_name:
          firstRangeLeft = -0.1
          firstRangeRight = 0.1
          tailSign = -1
        elif "res_cotth" in res_hist_name:
          firstRangeLeft = -0.01
          firstRangeRight = 0.01
          tailSign = -1
        elif "res_phi" in res_hist_name:
          firstRangeLeft = -0.005
          firstRangeRight = 0.005
          tailSign = -1
        elif "res_pt" in res_hist_name:
          firstRangeLeft = -1.
          firstRangeRight = 2.0
          tailSign = -1
          isPtPlot = True
        else:
          print "Initial fit parameters not specified -- check resolution hist names!"
          firstRangeLeft = -2.
          firstRangeRight = 2.0
          tailSign = 0

        temp_gaus = ROOT.TF1("temp_gaus","gaus", firstRangeLeft, firstRangeRight )

        if do_control_fit_plots:
          c = ROOT.TCanvas("c","c", 600, 600)

        r = temp_res.Fit(temp_gaus, "SRML Q")

        if do_control_fit_plots:
#          r.Draw() #save fits for checks
          c.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/13TeV_011014_res/fit_checks_step1/fit" + res_hist_name + "_bin" + str(i) + ".png")
          c.Close()


        tmp_mean = r.Parameter(1)
        tmp_sigma = r.Parameter(2)
        xmin = tmp_mean - tmp_sigma*5 # update mean and sigma according to the first fit
        xmax = tmp_mean + tmp_sigma*5

        x = ROOT.RooRealVar("x","x", xmin, xmax)


        if tmp_mean < 0:
          meanRangeMin = 5*tmp_mean
          meanRangeMax = -5*tmp_mean
        else:
          meanRangeMin = -5*tmp_mean
          meanRangeMax = 5*tmp_mean

        meanRoo = ROOT.RooRealVar("mean","mean of gaussian", tmp_mean, meanRangeMin, meanRangeMax)
        sigmaRoo = ROOT.RooRealVar("sigma", "width of gaussian", tmp_sigma, tmp_sigma*0.5, tmp_sigma*1.5)


        if isPtPlot:
            meanRoo.setRange(firstRangeLeft, firstRangeRight)
            x.setRange(firstRangeLeft, firstRangeRight)

            sigmaRoo.setMin(tmp_sigma*0.1)
            sigmaRoo.setMax(tmp_sigma*5.)
        
        # CB parameters
        a = ROOT.RooRealVar("a","a",3.,0.3,10.)
        aDx = ROOT.RooRealVar("aDx","aDx",3.,1.,10.)
        n = ROOT.RooRealVar("n","n",5.,0.,10.)
        nDx = ROOT.RooRealVar("nDx","nDx",5.,0.,10.)


#        if tailSign != 0:
#            if tailSign < 0:
#                a.setVal(1); a.setMin(0.1); a.setMax(3.)
#                x.setRange(firstRangeLeft,firstRangeRight/3.) 
#            else:
#                aDx.setVal(1); aDx.setMin(0.1); aDx.setMax(3.)
#                x.setRange(firstRangeLeft/3.,firstRangeRight) 
            
        if i>1:
            meanRoo.setVal(previousMean)
            sigmaRoo.setVal(previousSigma)
            a.setVal(previousA)
            n.setVal(previousN)
            aDx.setVal(previousADx)
            nDx.setVal(previousNDx)
            
        f_cb = ROOT.RooDoubleCB("cb","cb PDF",x,meanRoo,sigmaRoo,a,n,aDx,nDx)

#       temp_res_sub = ROOT.TH1Subset(temp_res, xmin, xmax)
        dh = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(x), ROOT.RooFit.Import(temp_res))

        if do_control_fit_plots:
            c2 = ROOT.TCanvas("c","c", 600, 600)
            xframe = x.frame(ROOT.RooFit.Title("Gaussian p.d.f."))
            
        f_cb.fitTo(dh, ROOT.RooFit.NumCPU(8), ROOT.RooFit.PrintLevel(-3)) #FIT HERE
        if do_control_fit_plots:
        
            dh.plotOn(xframe)
            f_cb.plotOn(xframe)
            xframe.Draw()

            c2.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/13TeV_011014_res/fit_checks_step2/fit" + res_hist_name + "_bin" + str(i) + ".png")
            c2.Close()


        tmp_mean = meanRoo.getVal()
        tmp_sigma = sigmaRoo.getVal()

        sigmaErr = sigmaRoo.getError()
        meanErr = meanRoo.getError()

        previousMean = meanRoo.getVal();
        previousSigma = sigmaRoo.getVal();
        previousA = a.getVal();
        previousN = n.getVal();
        previousADx = aDx.getVal();
        previousNDx = nDx.getVal();




        #-------- get 68% and 95% ------------
        fullAverage = temp_res.GetMean()
        step = temp_res.GetBinWidth(1)
        peakBin = temp_res.FindBin(tmp_mean)
        fullIntegral = temp_res.Integral(0, temp_res.GetNbinsX()+1)
        fullIntegralSx = temp_res.Integral(0, peakBin)
        fullIntegralDx = temp_res.Integral(peakBin, temp_res.GetNbinsX()+1)

        found68=False
        found95=False

        range68=0
        range95=0

        uncert68=0
        uncert95=0

        print "fullIntegral = " + str(fullIntegral)
        for j in range(0,temp_res.GetNbinsX()/2+1):
            fraction=temp_res.Integral(peakBin-j, peakBin+j)/fullIntegral
            fractionSx=temp_res.Integral(peakBin-j, peakBin)/fullIntegralSx
            fractionDx=temp_res.Integral(peakBin, peakBin+j)/fullIntegralDx

            if found68 and found95:
                break

            averageBinContent = (temp_res.GetBinContent(peakBin-j)+temp_res.GetBinContent(peakBin+j))/2
#            print "fraction = " + str(fraction) ## for debug

            if fraction > 0.68 and not found68:
                print "pass 0.68"
                found68=True            
                range68 = step*(2*j+1)*0.5

                print "average bin content = " + str(averageBinContent)
                print "average step = " + str(step)
                print "fullIntegral = " + str(fullIntegral)
                try:
                    uncert68 = math.sqrt(0.682*(1-0.682)/fullIntegral)/(averageBinContent/step/fullIntegral)
                except ZeroDivisionError:
                    uncert68 = 0.00001

            if fraction > 0.95 and not found95:    
                print "pass 0.95"
                found95=True
                range95 = step*(2*j+1)*0.5
                try:
                    uncert95 = math.sqrt(0.9*(1-0.9)/fullIntegral)/(averageBinContent/step/fullIntegral)
                except ZeroDivisionError:
                    uncert95 = 0.00001

        print "res 95 = " + str(range95) + ",  res68 = " + str(range68) 
        if mode == "68":
            res_1d.SetBinContent(i, range68)
            res_1d.SetBinError(i, uncert68)

        if mode == "95":
            res_1d.SetBinContent(i, range95)
            res_1d.SetBinError(i, uncert95)

    
    return res_1d    
