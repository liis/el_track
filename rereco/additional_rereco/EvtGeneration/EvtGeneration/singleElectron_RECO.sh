

#cmsDriver.py  --conditions auto:run2_mc -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT,RAW2DIGI,L1Reco -n 10 --magField 38T_PostLS1 --eventcontent FEVTDEBUGHLT --io ./outfiles/SingleElectronPt10_RECO.io --python SingleElectronPt10_RECO.py --no_exec --filein file:./outfiles/SingleElectronPt10_GEN.root --fileout outfiles/file:SingleElectronPt10_RECO.root 

#cmsDriver.py --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --conditions auto:run2_mc -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@relval,RAW2DIGI,L1Reco --datatier GEN-SIM-DIGI-RAW-HLTDEBUG -n 10 --magField 38T_PostLS1 --eventcontent FEVTDEBUGHLT --io DIGIUP15.io --python DIGIUP15.py --no_exec --filein file:./outfiles/SingleElectronPt10_GEN.root --fileout file:step2.root

#cmsDriver.py --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --conditions auto:run2_mc -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT,RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-RECO-RECODEBUG -n 10 --magField 38T_PostLS1 --eventcontent RECODEBUG --io DIGIUP15.io --python DIGIUP15.py --no_exec --filein file:./outfiles/SingleElectronPt10_GEN.root --fileout file:step2.root

#-------- testing run from RAW -----------------
cmsDriver.py --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --conditions PRE_STA71_V4::All -s RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-RECO-RECODEBUG -n -1 --magField 38T_PostLS1 --eventcontent RECODEBUG --io ./outfiles/testing/SingleElectronPt10_RECO.io --python ./conf_files_cmsdr/testing/SingleElectronPt10_RECO.py --filein file:./outfiles/testing/SingleElectronPt10_RAW.root --fileout file:./outfiles/testing/SingleElectronPt10_RECO.root 