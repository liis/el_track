

cmsDriver.py SingleElectronPt10_cfi --conditions auto:run2_mc -n 10 --eventcontent FEVTDEBUG --relval 9000,500 -s GEN,SIM --datatier GEN-SIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --magField 38T_PostLS1 --io ./outfiles/SingleElectronPt10_GEN.io --python ./conf_files_cmsdr/SingleElectronPt10_GEN.py --no_exec --fileout file:./outfiles/SingleElectronPt10_GEN.root

#cmsDriver.py SingleElectronPt10_cfi --conditions auto:run2_mc -n 10 --eventcontent RECODEBUG -s GEN,SIM --magField 38T_PostLS1 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --io ./outfiles/SingleElectronPt10_GEN.io --python SingleElectronPt10_GEN.py --no_exec --fileout file:./outfiles/SingleElectronPt10_GEN.root