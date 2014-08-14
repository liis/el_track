# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --conditions auto:run2_mc -s HLT,RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-RECO-RECODEBUG -n 10 --magField 38T_PostLS1 --eventcontent RECODEBUG --io DIGIUP15.io --python DIGIUP15.py --no_exec --filein /store/user/liis/GSF_tracking_samples/test_Zee_62_RECO_2/044C97C0-B675-E311-B49C-0025901D4D76.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/user/liis/GSF_tracking_samples/test_Zee_62_RECO_2/044C97C0-B675-E311-B49C-0025901D4D76.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('--customise nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECODEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:step2.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO-RECODEBUG')
    )
)

from SimGeneral.MixingModule.aliases_cfi import *
from SimGeneral.MixingModule.mixObjects_cfi import *
from SimGeneral.MixingModule.trackingTruthProducer_cfi import *


process.mix = cms.EDProducer("MixingModule",
                             digitizers = cms.PSet(
    mergedtruth = cms.PSet(
    trackingParticles
    )
    ),
                             LabelPlayback = cms.string(' '),
                             maxBunch = cms.int32(3),
                             minBunch = cms.int32(-5), ## in terms of 25 ns
                             
                             bunchspace = cms.int32(25),
                             mixProdStep1 = cms.bool(False),
                             mixProdStep2 = cms.bool(False),
                                 
                             playback = cms.untracked.bool(False),
                             useCurrentProcessOnly = cms.bool(False),
                             mixObjects = cms.PSet(
    mixTracks = cms.PSet(
    mixSimTracks
    ),
    mixVertices = cms.PSet(
    mixSimVertices
    ),
    mixSH = cms.PSet(
    mixSimHits
    ),
    mixHepMC = cms.PSet(
    mixHepMCProducts
    )
    )
                             )


process.p1 = cms.Path(process.mix)
# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule()
#process.schedule.extend(process.p1)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.p1, process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECODEBUGoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
