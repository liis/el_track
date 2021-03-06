from array import *
import numpy as np
maxpart = 300
maxhit = 30

int_list = ["np_reco", "np_gen", "np_fake", "np_gen_toReco"]
int_array_list = ["gen_nrMatchedSeedHits", "gen_nrSimHits", "gen_nrUniqueSimHits", "gen_nrMatchedRecHits","gen_reco_matched", "gen_reco_matched_sim", "charge_matched", "is_ecalDrivenSeed", "is_trackerDrivenSeed", "gen_matchedSeedOkCharge"] 

int_2D_array_list = ["gen_hit_subdetector", "gen_hit_layer"]
float_2D_array_list = ["gen_hit_pt"]

double_list = []
double_array_list_hits = []
double_array_list_res_pull = ["pt_pull", "qoverp_pull", "theta_pull", "phi_pull", "dxy_pull", "dz_pull", "gen_matched_rec_pt", "gen_matched_cotth", "gen_matched_rec_cotth", "gen_matched_rec_qoverp", "gen_matched_qoverp", "gen_matched_dz", "gen_matched_rec_dz", "gen_matched_rec_phi", "gen_matched_dxy", "gen_matched_rec_dxy", "gen_matched_phi", "gen_mathed_rec_phi"] #, "gen_matched_theta", "gen_matched_rec_theta"]

eff_double_array_list = ["reco_pt", "reco_eta", "gen_pt", "gen_eta", "fake_pt", "fake_eta", "gen_matched_eta", "gen_matched_pt", "gen_matchedSeedQuality"]
eff_int_array_list = ["gen_reco_matched", "gen_reco_matched_sim", "charge_matched", "is_ecalDrivenSeed", "is_trackerDrivenSeed", "gen_nrMatchedSeedHits", "gen_matchedSeedOkCharge"]

double_array_list = ["reco_phi", "gen_phi", "fake_phi", "gen_matched_phi", "gen_dxy", "gen_dz", "gen_bremFraction"] + double_array_list_hits + double_array_list_res_pull + eff_double_array_list
int_array_list = ["gen_nrUniqueSimHits"] + eff_int_array_list 

def var_type(var_name):
    """
    return the variable type for tree entry
    """

    if var_name in int_list: return array('i',[0])
    if var_name in double_list: return array('d', [0])
    if var_name in int_array_list: return array('i',[0]*maxpart)
    if var_name in double_array_list: return array('d', [0]*maxpart)
    if var_name in int_2D_array_list: return np.array( np.zeros( (maxpart, maxhit) ), dtype=np.int32) 
    if var_name in float_2D_array_list: return np.array( np.zeros( (maxpart, maxhit) ), dtype=np.float32 ) 


eff_list = int_list + eff_double_array_list + eff_int_array_list 
var_list = int_list + int_array_list + double_list + double_array_list + int_2D_array_list + float_2D_array_list
