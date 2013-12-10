import ROOT
from array import array
from tree_variables import var_list, var_type, maxhit
from histlib import fill_hist_ratio, log_binning

debug = False

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ptmode', dest='ptmode', choices=["Pt1","Pt10","Pt100","FlatPt"], required=True, help= "pt cut in analyzed dataset")
parser.add_argument('--isGSF', dest='is_gsf',  action="store_true")
args = parser.parse_args()

indir = "$WORKING_DIR/tree_to_histo/input_trees/"

if args.is_gsf:
    infile = "trackValTree_el" + args.ptmode + "GSF.root"
else:
    infile = "trackValTree_el" + args.ptmode + ".root"

f = ROOT.TFile(indir + infile)
t = f.Get("TrackValTreeMaker/trackValTree")

report_every = 100
max_event = -1
nEvt = t.GetEntries()

print "Analyzing dataset for " + args.ptmode + " with " + str(nEvt) + " events"

#### initialize tree #######
vt = dict([ (v, var_type(v)) for v in var_list ]) #associate proper data type for variables in the tree

t.SetBranchStatus("*",0)
for v in var_list: # relate vList to the tree
    t.SetBranchStatus(v, 1)
    t.SetBranchAddress(v, vt[v])
    t.AddBranchToCache(v, ROOT.kTRUE)
    
t.StopCacheLearningPhase()

neta = 100
mineta = -2.5
maxeta = 2.5

npt = 50
minpt = 0.1
maxpt = 200

xbinspt = log_binning(npt,minpt,maxpt) #log binning for pt histograms  

h_reco_eta = ROOT.TH1F("reco_eta","nrReco vs eta", neta, mineta, maxeta)
h_fake_eta = ROOT.TH1F("fake_eta","nrFake vs eta", neta, mineta, maxeta)
h_sim_eta = ROOT.TH1F("sim_eta", "nrSim vs eta", neta, mineta, maxeta)
h_sim_to_reco_match_eta = ROOT.TH1F("nsimtoreco_eta", "nr sim to reco matched vs eta", neta,mineta,maxeta)

h_reco_pt_barrel = ROOT.TH1F("nreco_pt_barrel","nrReco vs pt barrel", npt, array('d',xbinspt) )
h_reco_pt_endcap = ROOT.TH1F("nreco_pt_endcap","nrReco vs pt endcap", npt, array('d',xbinspt) )
h_reco_pt_trans = ROOT.TH1F("nreco_pt_trans","nrReco vs pt trans", npt, array('d',xbinspt) )

h_fake_pt_barrel = ROOT.TH1F("nfake_pt_barrel","nrFake vs pt barrel", npt, array('d', xbinspt) )
h_fake_pt_endcap = ROOT.TH1F("nfake_pt_endcap","nrFake vs pt endcap", npt, array('d', xbinspt) )
h_fake_pt_trans = ROOT.TH1F("nfake_pt_trans","nrFake vs pt trans", npt, array('d',xbinspt) )

h_sim_pt_barrel = ROOT.TH1F("sim_pt_barrel", "nrSim vs pt", npt, array('d',xbinspt))
h_sim_pt_endcap = ROOT.TH1F("sim_pt_endcap", "nrSim vs pt", npt, array('d',xbinspt))
h_sim_pt_trans = ROOT.TH1F("sim_pt_trans", "nrSim vs pt", npt, array('d', xbinspt))

h_sim_to_reco_match_pt_barrel = ROOT.TH1F("nsimtoreco_pt_barrel", "nr sim to reco matched vs pt", npt,array('d', xbinspt))
h_sim_to_reco_match_pt_endcap = ROOT.TH1F("nsimtoreco_pt_endcap", "nr sim to reco matched vs pt", npt, array('d',xbinspt))
h_sim_to_reco_match_pt_trans = ROOT.TH1F("nsimtoreco_pt_trans", "nr sim to reco matched vs pt", npt, array('d', xbinspt))

sim_hit_quality_flags = {"quality065": 0.65,"quality075": 0.75,"quality085": 0.85, "quality095": 0.95}

h_sim_eta_byhit = []
h_sim_pt_byhit_barrel = []
h_sim_pt_byhit_trans = []
h_sim_pt_byhit_endcap = []

h_sim_eta_byhit_quality = {}
h_sim_pt_byhit_quality_barrel = {}
h_sim_pt_byhit_quality_trans = {}
h_sim_pt_byhit_quality_endcap = {}

#------------------------tracker-hit-histograms---------------------

for nrhit in range(0,maxhit):
    h_sim_eta_byhit.append(ROOT.TH1F("eta_at_"+str(nrhit), "eta_at_"+str(nrhit),neta,mineta,maxeta) )
    
    h_sim_pt_byhit_barrel.append(ROOT.TH1F("pt_at_"+str(nrhit)+"_barrel", "pt_at_"+str(nrhit),npt,array('d',xbinspt)) )
    h_sim_pt_byhit_trans.append(ROOT.TH1F("pt_at_"+str(nrhit)+"_trans", "pt_at_"+str(nrhit),npt, array('d', xbinspt)) )
    h_sim_pt_byhit_endcap.append(ROOT.TH1F("pt_at_"+str(nrhit)+"_endcap", "pt_at_"+str(nrhit),npt, array('d', xbinspt)) )



for quality_flag in sim_hit_quality_flags:
    h_sim_eta_byhit_quality[quality_flag] = []
    
    h_sim_pt_byhit_quality_barrel[quality_flag] = []
    h_sim_pt_byhit_quality_endcap[quality_flag] = []
    h_sim_pt_byhit_quality_trans[quality_flag] = []
    for nrhit in range(0,maxhit):
        h_sim_eta_byhit_quality[quality_flag].append(ROOT.TH1F("eta_at_"+str(nrhit)+"_"+quality_flag,"eta_at_"+str(nrhit)+"_"+quality_flag, neta, mineta, maxeta) ) #pt for tracks, that have nth track (save for efficiencies)
        h_sim_pt_byhit_quality_barrel[quality_flag].append(ROOT.TH1F("pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", "pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", npt, array('d', xbinspt)) )
        h_sim_pt_byhit_quality_endcap[quality_flag].append(ROOT.TH1F("pt_at"+str(nrhit)+"_"+quality_flag+"_endcap", "pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", npt, array('d', xbinspt)) )
        h_sim_pt_byhit_quality_trans[quality_flag].append(ROOT.TH1F("pt_at"+str(nrhit)+"_"+quality_flag+"_trans", "pt_at"+str(nrhit)+"_"+quality_flag+"_barrel", npt, array('d',xbinspt)) )
#--------------------------------------------------------------------


# Event loop
for i in range(nEvt):
    if i % report_every == 0: print "Event nr: " + str(i)
    if i == max_event and not i == -1: break

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

    for it_p in range( vt['np_gen'][0]): #loop over simulated tracks
        gen_track_pt = vt['gen_pt'][it_p]
        gen_track_eta = vt['gen_eta'][it_p]

        h_sim_eta.Fill(gen_track_eta)
        
#        print "Found generated track with pt " + str(gen_track_pt)
        for nrhit in range(0,maxhit):
            pt_at_entry = vt["gen_hit_pt"][it_p][nrhit]
            hit_subdet = vt["gen_hit_subdetector"][it_p][nrhit]
            hit_layer = vt["gen_hit_layer"][it_p][nrhit]

#            print "hit pt = " + str(pt_at_entry) + ", hit subdetector = " + str(hit_subdet) + ", layer = " + str(hit_layer);
            if pt_at_entry > 0: # Nth hit exists -- fill eta and pt histograms for efficiency plots
                h_sim_eta_byhit[nrhit].Fill(gen_track_eta)
                
                if abs(vt['gen_eta'][it_p]) < 0.9:
                    h_sim_pt_byhit_barrel[nrhit].Fill(gen_track_pt)
                elif abs(vt['gen_eta'][it_p]) < 1.4:
                    h_sim_pt_byhit_trans[nrhit].Fill(gen_track_pt)
                elif abs(vt['gen_eta'][it_p]) < 2.5:
                    h_sim_pt_byhit_endcap[nrhit].Fill(gen_track_pt)

                for quality_flag, cut in sim_hit_quality_flags.iteritems():
                    if pt_at_entry/gen_track_pt > cut: # check hit pt-quality condition at each hit 
                        h_sim_eta_byhit_quality[quality_flag][nrhit].Fill(gen_track_eta)

                        if abs(vt['gen_eta'][it_p]) < 0.9:
                            h_sim_pt_byhit_quality_barrel[quality_flag][nrhit].Fill(gen_track_pt) 
                        elif abs(vt['gen_eta'][it_p]) < 1.4:
                            h_sim_pt_byhit_quality_trans[quality_flag][nrhit].Fill(gen_track_pt) 
                        elif abs(vt['gen_eta'][it_p]) < 2.5:
                            h_sim_pt_byhit_quality_endcap[quality_flag][nrhit].Fill(gen_track_pt) 

                
        if abs(vt['gen_eta'][it_p]) < 0.9:
            h_sim_pt_barrel.Fill(vt['gen_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 1.4:
            h_sim_pt_trans.Fill(vt['gen_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 2.5:
            h_sim_pt_endcap.Fill(vt['gen_pt'][it_p])
            
    for it_p in range( vt['np_gen_toReco'][0]): # loop over matched simulated tracks
        h_sim_to_reco_match_eta.Fill(vt['gen_matched_eta'][it_p])

        if abs(vt['gen_eta'][it_p]) < 0.9:
            h_sim_to_reco_match_pt_barrel.Fill(vt['gen_matched_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 1.4:
            h_sim_to_reco_match_pt_trans.Fill(vt['gen_matched_pt'][it_p])
        elif abs(vt['gen_eta'][it_p]) < 2.5:
            h_sim_to_reco_match_pt_endcap.Fill(vt['gen_matched_pt'][it_p])

#----------------------------------efficiencies and fake rates wrt eta---------------------------------
h_fakerate_eta = fill_hist_ratio(h_fake_eta, h_reco_eta,"fake_rate_eta")
h_eff_eta = fill_hist_ratio(h_sim_to_reco_match_eta, h_sim_eta, "eff_eta")    

h_eff_hit_eta = {}
h_eff_hit_pt_barrel = {}
h_eff_hit_pt_endcap = {}
h_eff_hit_pt_trans = {}
for quality_flag in sim_hit_quality_flags:
    h_eff_hit_eta[quality_flag] = []
    h_eff_hit_pt_barrel[quality_flag] = []
    h_eff_hit_pt_trans[quality_flag] = []
    h_eff_hit_pt_endcap[quality_flag] = []

    for nrhit in range(0,maxhit):
        h_eff_hit_eta[quality_flag].append( fill_hist_ratio(h_sim_eta_byhit_quality[quality_flag][nrhit], h_sim_eta_byhit[nrhit], "eff_eta_athit_"+str(nrhit)+"_"+quality_flag) )
        
        h_eff_hit_pt_barrel[quality_flag].append( fill_hist_ratio(h_sim_pt_byhit_quality_barrel[quality_flag][nrhit], h_sim_pt_byhit_barrel[nrhit], "eff_pt_athit_"+str(nrhit)+"_barrel_"+quality_flag, binning = "log") )
        h_eff_hit_pt_trans[quality_flag].append( fill_hist_ratio(h_sim_pt_byhit_quality_trans[quality_flag][nrhit], h_sim_pt_byhit_trans[nrhit], "eff_pt_athit_"+str(nrhit)+"_trans_"+quality_flag, binning = "log") )
        h_eff_hit_pt_endcap[quality_flag].append( fill_hist_ratio(h_sim_pt_byhit_quality_endcap[quality_flag][nrhit], h_sim_pt_byhit_endcap[nrhit], "eff_pt_athit_"+str(nrhit)+"_endcap_"+quality_flag, binning = "log") )

# efficiency and fake rate wrt pt
h_eff_pt_barrel = fill_hist_ratio(h_sim_to_reco_match_pt_barrel, h_sim_pt_barrel, "eff_pt_barrel", binning="log")
h_eff_pt_trans = fill_hist_ratio(h_sim_to_reco_match_pt_trans, h_sim_pt_trans, "eff_pt_trans", binning="log")
h_eff_pt_endcap = fill_hist_ratio(h_sim_to_reco_match_pt_endcap, h_sim_pt_endcap, "eff_pt_endcap", binning="log")

h_fakerate_pt_barrel = fill_hist_ratio(h_fake_pt_barrel, h_reco_pt_barrel, "fake_rate_pt_barrel", binning = "log")
h_fakerate_pt_trans = fill_hist_ratio(h_fake_pt_trans, h_reco_pt_trans, "fake_rate_pt_trans", binning = "log")
h_fakerate_pt_endcap = fill_hist_ratio(h_fake_pt_endcap, h_reco_pt_endcap, "fake_rate_pt_endcap", binning = "log")

if args.is_gsf:
    outfile = "histograms/trackValHistograms" + args.ptmode + "GSF.root"
else:
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

#-----efficiencies-----------
h_fakerate_eta.Write()
h_eff_eta.Write()

h_fakerate_pt_barrel.Write()
h_fakerate_pt_endcap.Write()
h_fakerate_pt_trans.Write()

h_eff_pt_barrel.Write()
h_eff_pt_endcap.Write()
h_eff_pt_trans.Write()

#------------hit-by-hit comparison------------
dir = o.mkdir("VariablesBySimhit")
dir.cd()
for quality_flag in sim_hit_quality_flags:
    for nrhit in range(0,maxhit):
        h_eff_hit_eta[quality_flag][nrhit].Write()
        h_sim_eta_byhit_quality[quality_flag][nrhit].Write()

        h_eff_hit_pt_barrel[quality_flag][nrhit].Write()
        h_eff_hit_pt_endcap[quality_flag][nrhit].Write()
        h_eff_hit_pt_trans[quality_flag][nrhit].Write()

        h_sim_pt_byhit_quality_barrel[quality_flag][nrhit].Write()
        h_sim_pt_byhit_quality_endcap[quality_flag][nrhit].Write()
        h_sim_pt_byhit_quality_trans[quality_flag][nrhit].Write()
for nrhit in range(0,maxhit):
    h_sim_eta_byhit[nrhit].Write()

    h_sim_pt_byhit_barrel[nrhit].Write()
    h_sim_pt_byhit_endcap[nrhit].Write()
    h_sim_pt_byhit_trans[nrhit].Write()

o.Close()
