
#--- real thing --
#cmsDriver.py SingleElectronPt10_cfi --conditions PRE_STA71_V4::All -n 100000 --eventcontent FEVTDEBUG --relval 9000,500 -s GEN,SIM --datatier GEN-SIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --magField 38T_PostLS1 --io ./outfiles/SingleElectronPt10_GEN.io --python ./conf_files_cmsdr/SingleEl_GEN/SingleElectronPt10_GEN.py --no_exec --fileout file:./outfiles/SingleElectronPt10_GEN.root

# --- testing ---
#cmsDriver.py SingleElectronPt100_cfi --conditions PRE_STA71_V4::All -n 100 --eventcontent FEVTDEBUG --relval 9000,500 -s GEN,SIM --datatier GEN-SIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --magField 38T_PostLS1 --io ./outfiles/testing/SingleElectronPt100_GEN.io --python ./conf_files_cmsdr/SingleElectronPt10_GEN.py --no_exec --fileout file:./outfiles/testing/SingleElectronPt100_GEN.root

# ------- test gen+raw --------
cmsDriver.py SingleElectronPt10_cfi --conditions PRE_STA71_V4::All -n 30 --eventcontent FEVTDEBUG --relval 9000,500 -s GEN,SIM,DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@relval,RAW2DIGI,L1Reco --datatier GEN-SIM-FEVTDEBUG --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --magField 38T_PostLS1 --io ./outfiles/testing/SingleElectronPt10_RAW.io --python ./conf_files_cmsdr/testing/SingleElectronPt10_RAW.py --fileout file:./outfiles/testing/SingleElectronPt10_RAW.root

#cmsDriver.py SingleElectronPt10_cfi --conditions auto:run2_mc -n 10 --eventcontent RECODEBUG -s GEN,SIM --magField 38T_PostLS1 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --io ./outfiles/SingleElectronPt10_GEN.io --python SingleElectronPt10_GEN.py --no_exec --fileout file:./outfiles/SingleElectronPt10_GEN.root