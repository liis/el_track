import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

process.source = cms.Source("PoolSource",
                            #         fileNames = cms.untracked.vstring(filePrefex+sys.argv[2]),
                            #     inputCommands = cms.untracked.vstring("drop *","keep *_source_*_*"),
                            fileNames = cms.untracked.vstring(
#        '/store/user/liis/GSF_tracking_samples/test_Zee_62_RECO_2/044C97C0-B675-E311-B49C-0025901D4D76.root'), #t3 ch
        'file:/tmp/liis/044C97C0-B675-E311-B49C-0025901D4D76.root' #lxplus0157
        ),
                            )

# Production Info
process.configurationMetadata = cms.untracked.PSet( # ??
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('step3 nevts:50'),
    name = cms.untracked.string('Applications')
    )


#process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
#                                           splitLevel = cms.untracked.int32(0),
#                                           SelectEvents = cms.untracked.PSet(  SelectEvents = ( cms.vstring( 'reconstruction_step',))), 
#                                           eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#                                           outputCommands = process.RECODEBUGEventContent.outputCommands,
#                                           fileName = cms.untracked.string('file:samtest_reco.root'),
#                                           dataset = cms.untracked.PSet(
#    filterName = cms.untracked.string(''),
#    dataTier = cms.untracked.string('GEN-SIM-RECO-RECODEBUG')
#    )
#)

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
                                         splitLevel = cms.untracked.int32(0),
                                         SelectEvents = cms.untracked.PSet(  SelectEvents = ( cms.vstring( 'reconstruction_step',))),
                                         eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
                                         outputCommands = process.RECODEBUGEventContent.outputCommands,
                                         fileName = cms.untracked.string('samtest_reco.root'),
                                         dataset = cms.untracked.PSet(
    filterName = cms.untracked.string(''),
    dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
                                         )
process.RECOSIMoutput.outputCommands.append('keep *PSimHits_g4SimHits_*_*')

# Other statements
process.mix.playback = True
#process.mix.digitizers = cms.PSet()
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
#for a in process.aliases: delattr(process, a) #?? #THIS IS GUILTY for the removal of "pixelDigiSimLink"

process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from CondCore.DBCommon.CondDBSetup_cfi import *
process.ecalES1 = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string("frontier://FrontierProd/CMS_COND_31X_ECAL"),
                               toGet = cms.VPSet(
    cms.PSet(record = cms.string("EcalChannelStatusRcd"),tag=cms.string("EcalChannelStatus_coll12_v1_mc")),
    cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),tag=cms.string("EcalIntercalibConstants_2011_V3_Bon_start_mc")),
    cms.PSet(record = cms.string("EcalLaserAlphasRcd"),tag=cms.string("EcalLaserAlphas_mc")),
    cms.PSet(record = cms.string("EcalPedestalsRcd"),tag=cms.string("EcalPedestals_v02_mc")),
    cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),tag=cms.string("EcalLaserAPDPNRatios_p1p2p3_v2_mc")),
    )
                               )

process.ecalES2 = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string("frontier://FrontierProd/CMS_COND_34X_ECAL"),
                               toGet = cms.VPSet(
    #  cms.PSet(record = cms.string("EcalSRSettingsRcd"),tag=cms.string("null")),
    cms.PSet(record = cms.string("EcalTPGLinearizationConstRcd"),tag=cms.string("EcalTPGLinearizationConst_beamv5_startup_mc")),
    )
                               )


process.ecalES4 = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string("frontier://FrontierProd/CMS_COND_ECAL_000"),
                               toGet = cms.VPSet(
    cms.PSet(record = cms.string("EcalLinearCorrectionsRcd"),tag=cms.string("EcalLinearCorrections_mc")),
    )
                               )


process.ecalES5 = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string("frontier://FrontierProd/CMS_COND_31X_PRESHOWER"),
                               toGet = cms.VPSet(
    cms.PSet(record = cms.string("ESChannelStatusRcd"),tag=cms.string("ESChannelStatus_LG_V04_mc")),
    )
                               )

process.es_prefer_ecal1 = cms.ESPrefer("PoolDBESSource","ecalES1")
process.es_prefer_ecal2 = cms.ESPrefer("PoolDBESSource","ecalES2")
#process.es_prefer_ecal3 = cms.ESPrefer("PoolDBESSource","ecalES3")
process.es_prefer_ecal4 = cms.ESPrefer("PoolDBESSource","ecalES4")
process.es_prefer_ecal5 = cms.ESPrefer("PoolDBESSource","ecalES5")

#process.mcFilter = cms.EDFilter("MCTruthFilter",
#                                genParticlesTag = cms.InputTag("genParticles"),
#                                pid=cms.int32(11)
#                                )

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
#process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)

process.mixing = cms.Path(process.mix)

# Schedule definition
process.schedule = cms.Schedule(process.mixing,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.endjob_step,process.RECOSIMoutput_step)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

print process.particleFlowClusterECAL.inputECAL

# End of customisation functions
#if process.particleFlowClusterECAL.inputECAL.getModuleLabel()=="particleFlowClusterECALWithTimeSelected":
#    print "3D Timing: ON 3DTiming"
#else:
 #   print "3D Timing: OFF"
