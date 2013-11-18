#voms-proxy-init -voms cms:/cms/Role=priorityuser
[CRAB]
jobtype                    = cmssw
scheduler                  = glite
use_server                 = 0
                          
[CMSSW]                   
first_lumi                 = 1
#datasetpath                = /SingleElMinusFlatLogPt_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-8de1ffbdb519f9edbafc5606a1926f13/USER
#datasetpath                =/SingleElMinusPt1_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-7bc796286602c18e9ed77a7f93a692b8/USER
#datasetpath                =/SingleElMinusPt10_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG-v2/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-7bc796286602c18e9ed77a7f93a692b8/USER
datasetpath                =/SingleElMinusPt100_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG-v2/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-7bc796286602c18e9ed77a7f93a692b8/USER

no_block_boundary = 1

pset                       = ./reTracking.py
events_per_job             = 5000
total_number_of_events     = 100000
output_file                = step2.root
dbs_url                    = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet
                          
[USER]                    
#ui_working_dir             = ./crab_0_SingleElMinusFlatPt_START-v4
#ui_working_dir             = ./crab_0_SingleElMinusPt1_START-v4
#ui_working_dir             = ./crab_0_SingleElMinusPt10_START-v4
ui_working_dir             = ./crab_0_SingleElMinusPt100_START-v4
return_data                = 0
copy_data                  = 1
storage_element            = T2_US_UCSD
local_stage_out            = 0
publish_data               = 1
#publish_data_name          = SingleElMinusFlatPt_CMSSW_4_2_8-START42_V12_reTrackingGsf-v4
#publish_data_name          = SingleElMinusPt1_CMSSW_4_2_8-START42_V12_reTrackingGsf-v4
#publish_data_name          = SingleElMinusPt10_CMSSW_4_2_8-START42_V12_reTrackingGsf-v4
publish_data_name          = SingleElMinusPt100_CMSSW_4_2_8-START42_V12_reTrackingGsf-v4
dbs_url_for_publication    = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet


[GRID]
#role = priorityuser
#ce_black_list = T3_*
se_white_list = T2_US_UCSD
