import ROOT
from array import array
from tree_variables import var_list, var_type

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ptmode', dest='ptmode', choices=["Pt1","Pt10","Pt100","FlatPt"], required=True, help= "pt cut in analyzed dataset")

args = parser.parse_args()

if args.ptmode == "Pt1":
    fixpt = 1
if args.ptmode == "Pt10":
    fixpt = 10
if args.ptmode == "Pt100":
    fixpt = 100
if args.ptmode == "FlatPt":
    fixpt = 1

indir = "$WORKING_DIR/tree_to_histo/input_trees/"
infile = "trackValTree_el" + args.ptmode + ".root"

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 10000
nEvt = t.GetEntries()

print "Analyzing dataset for " + args.ptmode + " with " + str(nEvt) + " events"

#### initialize tree #######
vt = dict([ (v, var_type(v)) for v in var_list ]) #associate proper data type for every variable in the tree

t.SetBranchStatus("*",0)
for v in var_list: # relate vList to the tree
    t.SetBranchStatus(v, 1)
    t.SetBranchAddress(v, vt[v])
    t.AddBranchToCache(v, ROOT.kTRUE)
    
t.StopCacheLearningPhase()

neta = 100
mineta = -2.5
maxeta = 2.5

h_reco_eta = ROOT.TH1F("nreco_eta","nrReco vs eta", neta, mineta, maxeta)
h_fake_eta = ROOT.TH1F("nfake_eta","nrFake vs eta", neta, mineta, maxeta)
h_sim_eta = ROOT.TH1F("nsim_eta", "nrSim vs eta", neta, mineta, maxeta)
h_sim_to_reco_match_eta = ROOT.TH1F("nsimtoreco_eta", "nr sim to reco matched vs eta", neta,mineta,maxeta)

h_fakerate_eta = ROOT.TH1F("fake_rate_eta","Fake rate", neta, mineta, maxeta)
h_eff_eta = ROOT.TH1F("eff_eta","Efficiency", neta, mineta, maxeta)

npt = 50
minpt = 0.1
maxpt = 200

logxmin = ROOT.TMath.log10(minpt); #get logarithmic bins for pt histograms
print "logxmin = " + str(logxmin)
logxmax = ROOT.TMath.log10(maxpt);
binwidth = (logxmax-logxmin)/npt;
xbins = []
xbins.append(minpt);
for i in range(1,npt+1):
    xbins.append(minpt + ROOT.TMath.Power(10,logxmin + i*binwidth))

#h_reco_pt_barrel = ROOT.TH1F("nreco_pt_barrel","nrReco vs pt barrel", npt, minpt, maxpt)
h_reco_pt_barrel = ROOT.TH1F("nreco_pt_barrel","nrReco vs pt barrel", npt, array('d',xbins) )
h_reco_pt_endcap = ROOT.TH1F("nreco_pt_endcap","nrReco vs pt endcap", npt, array('d',xbins) )
h_reco_pt_trans = ROOT.TH1F("nreco_pt_trans","nrReco vs pt trans", npt, array('d',xbins) )

#h_fake_pt_barrel = ROOT.TH1F("nfake_pt_barrel","nrFake vs pt barrel", npt, minpt, maxpt)
h_fake_pt_barrel = ROOT.TH1F("nfake_pt_barrel","nrFake vs pt barrel", npt, array('d', xbins) )
h_fake_pt_endcap = ROOT.TH1F("nfake_pt_endcap","nrFake vs pt endcap", npt, array('d', xbins) )
h_fake_pt_trans = ROOT.TH1F("nfake_pt_trans","nrFake vs pt trans", npt, array('d',xbins) )

h_sim_pt_barrel = ROOT.TH1F("sim_pt_barrel", "nrSim vs pt", npt, array('d',xbins))
h_sim_pt_endcap = ROOT.TH1F("sim_pt_endcap", "nrSim vs pt", npt, array('d',xbins))
h_sim_pt_trans = ROOT.TH1F("sim_pt_trans", "nrSim vs pt", npt, array('d', xbins))

h_sim_to_reco_match_pt_barrel = ROOT.TH1F("nsimtoreco_pt_barrel", "nr sim to reco matched vs pt", npt,array('d', xbins))
h_sim_to_reco_match_pt_endcap = ROOT.TH1F("nsimtoreco_pt_endcap", "nr sim to reco matched vs pt", npt, array('d',xbins))
h_sim_to_reco_match_pt_trans = ROOT.TH1F("nsimtoreco_pt_trans", "nr sim to reco matched vs pt", npt, array('d', xbins))

#h_fakerate_pt_barrel = ROOT.TH1F("fake_rate_pt_barrel","Fake rate barrel", npt, minpt, maxpt)
h_fakerate_pt_barrel = ROOT.TH1F("fake_rate_pt_barrel","Fake rate barrel", npt, array('d',xbins))
h_fakerate_pt_endcap = ROOT.TH1F("fake_rate_pt_endcap","Fake rate endcap", npt, array('d',xbins))
h_fakerate_pt_trans = ROOT.TH1F("fake_rate_pt_trans","Fake rate trans", npt, array('d',xbins))

h_eff_pt_barrel = ROOT.TH1F("eff_pt_barrel","Efficiency barrel", npt, array('d',xbins))
h_eff_pt_endcap = ROOT.TH1F("eff_pt_endcap","Efficiency endcap", npt, array('d',xbins))
h_eff_pt_trans = ROOT.TH1F("eff_pt_trans","Efficiency trans", npt, array('d',xbins))

# Event loop
for i in range(nEvt):
    if i % report_every == 0: print "Event nr: " + str(i)
    t.LoadTree(i)
    t.GetEntry(i)
    for it_p in range( vt['np_reco'][0]): # loop over reco-tracks
        h_reco_eta.Fill(vt['reco_eta'][it_p])
        
        if abs(vt['reco_eta'][it_p]) < 0.9:
            h_reco_pt_barrel.Fill(vt['reco_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 1.4:
            h_reco_pt_trans.Fill(vt['reco_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 2.5:
            h_reco_pt_endcap.Fill(vt['reco_pt'][it_p])

    for it_p in range( vt['np_fake'][0]): # loop over fake traks
        h_fake_eta.Fill(vt['fake_eta'][it_p])

        if abs(vt['reco_eta'][it_p]) < 0.9:
            h_fake_pt_barrel.Fill(vt['fake_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 1.4:
            h_fake_pt_trans.Fill(vt['fake_pt'][it_p])
        elif abs(vt['reco_eta'][it_p]) < 2.5:
            h_fake_pt_endcap.Fill(vt['fake_pt'][it_p])

    for it_p in range( vt['np_gen'][0]):
        h_sim_eta.Fill(vt['gen_eta'][it_p])
        
        if abs(vt['gen_eta'][it_p]) < 0.9:
            h_sim_pt_barrel.Fill(vt['gen_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 1.4:
            h_sim_pt_trans.Fill(vt['gen_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 2.5:
            h_sim_pt_endcap.Fill(vt['gen_pt'][it_p])

    for it_p in range( vt['np_gen_toReco'][0]):
        h_sim_to_reco_match_eta.Fill(vt['gen_matched_eta'][it_p])

        if abs(vt['gen_eta'][it_p]) < 0.9:
            h_sim_to_reco_match_pt_barrel.Fill(vt['gen_matched_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 1.4:
            h_sim_to_reco_match_pt_trans.Fill(vt['gen_matched_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 2.5:
            h_sim_to_reco_match_pt_endcap.Fill(vt['gen_matched_pt'][it_p])

# efficiency and fake rate wrt eta
for i in range(1,neta+1): 
    if not h_reco_eta.GetBinContent(i) == 0:
        h_fakerate_eta.SetBinContent(i, h_fake_eta.GetBinContent(i)/h_reco_eta.GetBinContent(i) )
    else:
        h_fakerate_eta.SetBinContent(i,0)
        
    if not h_sim_eta.GetBinContent(i) == 0:
        h_eff_eta.SetBinContent(i, h_sim_to_reco_match_eta.GetBinContent(i)/h_sim_eta.GetBinContent(i) )
    else:
        h_sim_eta.SetBinContent(i,0)

# efficiency and fake rate wrt pt
for i in range(1,npt+1):
    if not h_reco_pt_barrel.GetBinContent(i) == 0:
        h_fakerate_pt_barrel.SetBinContent(i, h_fake_pt_barrel.GetBinContent(i)/h_reco_pt_barrel.GetBinContent(i) )
    else:
        h_fakerate_pt_barrel.SetBinContent(i,0)
    if not h_reco_pt_endcap.GetBinContent(i) == 0:
        h_fakerate_pt_endcap.SetBinContent(i, h_fake_pt_endcap.GetBinContent(i)/h_reco_pt_endcap.GetBinContent(i) )
    else:
        h_fakerate_pt_endcap.SetBinContent(i,0)

    if not h_reco_pt_trans.GetBinContent(i) ==0:
        h_fakerate_pt_trans.SetBinContent(i, h_fake_pt_trans.GetBinContent(i)/h_reco_pt_trans.GetBinContent(i) )
    else:
        h_fakerate_pt_trans.SetBinContent(i,0)
        
    if not h_sim_pt_barrel.GetBinContent(i) == 0:
        h_eff_pt_barrel.SetBinContent(i, h_sim_to_reco_match_pt_barrel.GetBinContent(i)/h_sim_pt_barrel.GetBinContent(i) )
    else:
        h_eff_pt_barrel.SetBinContent(i,0)
    if not h_sim_pt_endcap.GetBinContent(i) == 0:
        h_eff_pt_endcap.SetBinContent(i, h_sim_to_reco_match_pt_endcap.GetBinContent(i)/h_sim_pt_endcap.GetBinContent(i) )
    else:
        h_eff_pt_endcap.SetBinContent(i,0)
    if not h_sim_pt_trans.GetBinContent(i) == 0:
        h_eff_pt_trans.SetBinContent(i, h_sim_to_reco_match_pt_trans.GetBinContent(i)/h_sim_pt_trans.GetBinContent(i) )
    else:
        h_eff_pt_trans.SetBinContent(i,0)

outfile = "histograms/trackValHistograms" + args.ptmode  + ".root"    
o = ROOT.TFile(outfile,"recreate")
h_reco_eta.Write()
h_fake_eta.Write()
h_sim_eta.Write()
h_sim_to_reco_match_eta.Write()

h_reco_pt_barrel.Write()
h_reco_pt_endcap.Write()
h_reco_pt_trans.Write()

h_fake_pt_barrel.Write()
h_fake_pt_endcap.Write()
h_fake_pt_trans.Write()

h_sim_pt_barrel.Write()
h_sim_pt_endcap.Write()
h_sim_pt_trans.Write()

h_sim_to_reco_match_pt_barrel.Write()
h_sim_to_reco_match_pt_trans.Write()
h_sim_to_reco_match_pt_endcap.Write()

h_fakerate_eta.Write()
h_eff_eta.Write()

h_fakerate_pt_barrel.Write()
h_fakerate_pt_endcap.Write()
h_fakerate_pt_trans.Write()

h_eff_pt_barrel.Write()
h_eff_pt_endcap.Write()
h_eff_pt_trans.Write()
o.Close()
