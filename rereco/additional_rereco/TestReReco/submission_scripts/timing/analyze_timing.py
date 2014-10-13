import os
import ROOT

indir = os.getcwd() + "/input_timing"
summary_file_name = "final_timing_eltrk_summary.txt"
cat_dir = "timing_maxCand_1_maxChi2_2000_nSigma_3_gsf_1412948895"

print "Analyzing files from: " + indir


def get_el_trk_timing( infilename, norm_reference = 0.3 ):
    f = file(inpath)
    for line in f.readlines():
        if "SeedGeneratorFromRegionHitsEDProducer:pixelLessStepSeeds" in line: #get reference value (unaffected by electron tracking) to normalize timing
            reference = float(line.split(':')[-1].split('ms')[0])
            norm = reference/norm_reference
#            print "normalization factor = " + str( norm )

    f.seek(0)
    for line in f: # loop again to have the normalization reference
        if "CkfTrackCandidateMaker:electronCkfTrackCandidates" in line:
            el_timing = float(line.split(':')[-1].split('ms')[0])
            el_timing_norm = el_timing*norm
#            print "el timing = " + str(el_timing)
#            print "normalized el timing = " + str(el_timing_norm)

    return el_timing_norm



timing_hist = ROOT.TH1F("timing","timing", 10, 0, 10)

ibin = 0
for dir in os.listdir(indir): # loop over timing output directories
    if dir.split('_')[0] == "timing":
        print "Get input from: " + dir
        inpath = indir + "/" + dir + "/" + summary_file_name
        el_timing_norm = get_el_trk_timing(inpath)
        print "el timing norm = " + str(el_timing_norm)
        
        #------------ Write timing to histogram ---------------
        binlabel = dir.split('timing_')[1].split('_gsf')[0]
        timing_hist.SetBinContent(ibin, el_timing_norm)
        timing_hist.GetXaxis().SetBinLabel(ibin, binlabel)
        
        ibin+=1

outfile = "timing_summary.root"
o = ROOT.TFile(outfile,"recreate")
timing_hist.Write()
o.Close()
