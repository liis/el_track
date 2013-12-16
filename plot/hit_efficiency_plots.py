import sys, os
import tdrstyle
tdrstyle.tdrstyle()
import ROOT
from plot_fun import infilenames_eta, infilenames_pt, load_input_files, get_hit_efficiency_hist, draw_efficiency_histograms, draw_legend

indir = "$WORKING_DIR/tree_to_histo/histograms/"

infiles_eta = load_input_files(indir, infilenames_eta)
infiles_pt = load_input_files(indir, infilenames_pt)

hist_hit_eff = {}
hist_hit_eff[ "eta_065"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 65, nrhits = 14)
hist_hit_eff["eta_075"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 75, nrhits = 14)
hist_hit_eff["eta_085"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 85, nrhits = 14)
hist_hit_eff["eta_095"] = get_hit_efficiency_hist(infiles_eta, "eta", quality = 95, nrhits = 14)

hist_hit_eff["pt_065"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 65, nrhits = 14)
hist_hit_eff["pt_075"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 75, nrhits = 14)
hist_hit_eff["pt_085"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 85, nrhits = 14)
hist_hit_eff["pt_095"] = get_hit_efficiency_hist(infiles_pt, "pt", quality = 95, nrhits = 14)


quality_points = [65,75,85,95]
eta_regions = ["barrel", "trans", "endcap"]

hit_hists = {}
for eta_region in eta_regions:
    hit_hists[eta_region] = {}
    for quality_point in quality_points:
        hist_name = "VariablesBySimhit/nr_hit_quality0" + str(quality_point) + "_" + eta_region
        print "loading histogram " + hist_name
        hit_hists[eta_region][quality_point] = infiles_pt["FlatPt"].Get(hist_name) 
        
print hit_hists["barrel"]

for eta_region, hists_by_quality in hit_hists.iteritems():
    c_eff = ROOT.TCanvas("c_eff_" + eta_region, eta_region)

    print hists_by_quality.values()
    print hists_by_quality.keys()
    draw_efficiency_histograms(hists_by_quality.values(), xtitle = "hit nr. in tracker (" + eta_region + ")")

    leg = ROOT.TLegend(0.6,0.7,0.9,0.9);
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    for quality in hists_by_quality.keys():
        leg.AddEntry(hists_by_quality[quality], "eff = " + str(quality) + "%")
    leg.Draw()

    c_eff.SaveAs("$WORKING_DIR/plot/out_plots/eff_by_hit_" + eta_region + ".pdf")
    c_eff.SaveAs("$WORKING_DIR/plot/out_plots/eff_by_hit_" + eta_region + ".png")


#---------------plots--------------------
for quality, hist_dictionary in hist_hit_eff.iteritems():
    for region, hists in hist_dictionary.iteritems():
        print "drawing histograms in region " + region

        c_eff = ROOT.TCanvas("c_eff_" + region + "_" +quality + "_" + region, "")
        c_eff.SetGrid()
        c_eff.SetRightMargin(0.25) #legend out of the plot area

        if region[:6] == "FlatPt":
            c_eff.SetLogx()

        draw_efficiency_histograms(hists,region = region, xtitle = "(quality = " + quality + ", pt = " + region + ")")
       
        leg = ROOT.TLegend(0.76,0.27,0.93,0.90); #out of the box coordinates
        leg.SetBorderSize(0)
        leg.SetFillColor(0)

        for ihit in range(0, len(hists) ):
            leg.AddEntry(hists[ihit], "tr. hit " + str(ihit + 1) )
        leg.Draw()

        c_eff.SaveAs("$WORKING_DIR/plot/out_plots/hiteff_" + quality + "_" + region + ".pdf")
        c_eff.SaveAs("$WORKING_DIR/plot/out_plots/hiteff_" + quality + "_" + region + ".png")
