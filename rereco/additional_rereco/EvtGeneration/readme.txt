Follow the sample configuration files of the existing singleElectron sample:/RelValSingleElectronPt35_UP15/CMSSW_7_1_0-POSTLS171_V15-v1/GEN-SIM-RECO

Get step1 cfg files from cmssw/Configuration/Generator/python, read automatically from this folder

Get cfg files from DAS:
Step 1. cmsDriver.py SingleElectronPt35_cfi --conditions auto:run2_mc -n 10 --eventcontent FEVTDEBUG --relval 9000,500 -s GEN,SIM --datatier GEN-SIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --magField 38T_PostLS1 --io SingleElectronPt35_UP15.io --python SingleElectronPt35_UP15.py --no_exec --fileout file:step1.root

Step 2: cmsDriver.py --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --conditions auto:run2_mc -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@relval,RAW2DIGI,L1Reco --datatier GEN-SIM-DIGI-RAW-HLTDEBUG -n 10 --magField 38T_PostLS1 --eventcontent FEVTDEBUGHLT --io DIGIUP15.io --python DIGIUP15.py --no_exec --filein file:step1.root --fileout file:step2.root

Step 3: cmsDriver.py --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --conditions auto:run2_mc -s RAW2DIGI,L1Reco,RECO,EI,VALIDATION,DQM --datatier GEN-SIM-RECO,DQM -n 10 --magField 38T_PostLS1 --eventcontent RECOSIM,DQM --io RECOUP15.io --python RECOUP15.py --no_exec --filein file:step2.root --fileout file:step3.root

Step 4: cmsDriver.py --conditions auto:run2_mc -s HARVESTING:validationHarvesting+dqmHarvesting --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --mc --magField 38T_PostLS1 --io HARVESTUP15.io --python HARVESTUP15.py -n 100 --no_exec --filein file:step3_inDQM.root --fileout file:step4.root


Magnetic field conditions from: /Configuration/StandardSequences/python/MagneticField_*
