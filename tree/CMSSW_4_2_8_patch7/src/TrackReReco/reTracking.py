import FWCore.ParameterSet.Config as cms

process = cms.Process("reGsfTracking")

# message logger
process.MessageLogger = cms.Service("MessageLogger",
     default = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

# source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
#readFiles.extend( ['file:/hdfs/cms/store/user/liis/El_GSF_studies/Pt100/step2_10_1_dDh.root'])
readFiles.extend( ['file:/hdfs/cms/store/user/liis/El_GSF_studies/test/SingleElMinusPt10_1_1_9Rd.root'])

process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START42_V11::All'

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")


maxCandDefault = 5
maxChi2Default = 2000
nSigmaDefault = 3.0

maxCand = 6
maxChi2 = 100
nSigma = 4

########################################################################
# to change parameters  as in slides from A.Tropiano

process.TrajectoryBuilderForElectrons.maxCand = cms.int32( maxCand )
process.ElectronChi2.MaxChi2 = cms.double( maxChi2 ) 
process.ElectronChi2.nSigma = cms.double( nSigma )

########################################################################

from RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_cfi import *
from RecoPixelVertexing.Configuration.RecoPixelVertexing_cff import *
from RecoEcal.Configuration.RecoEcal_cff import *
from RecoLocalCalo.Configuration.hcalGlobalReco_cff import *


process.load("SimG4Core.Application.g4SimHits_cfi")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

# Reconstruction geometry service
#process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
# geometry (Only Ecal)
#process.load("Geometry.EcalCommonData.EcalOnly_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

# use trivial ESProducer for tests
#process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")

# ECAL 'slow' digitization sequence
process.load("SimCalorimetry.Configuration.ecalDigiSequence_cff")

#process.load("SimCalorimetry.Configuration.ecalDigiSequenceComplete_cff")
#from SimCalorimetry.Configuration.ecalDigiSequence_cff import ecalDigiSequence


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   moduleSeeds = cms.PSet(
    g4SimHits = cms.untracked.uint32(9876),
    simEcalUnsuppressedDigis = cms.untracked.uint32(12345)
    )
                                                   )


# sequence for re-running gsfTracking over RECO
process.myGsfReco = cms.Sequence(
    process.offlineBeamSpot+
    process.siPixelDigis*process.siPixelClusters*
    process.siStripDigis*process.siStripZeroSuppression*process.siStripClusters*
    process.siPixelRecHits*
    process.recopixelvertexing*
   
    process.siStripMatchedRecHits* #make local hits
    process.iterTracking*process.trackCollectionMerging*  #CTF iterative tracking
    process.newCombinedSeeds* #together with EG Clasters, input for ecalDriven seeds

    process.g4SimHits*
    process.mix*
    process.ecalDigis*
    process.ecalPreshowerDigis*
#    process.ecalDigiSequence
    process.ecalLocalRecoSequence* # contains process.ecalRecHit

#    process.hbheprereco*
#    process.hcalGlobalRecoSequence*
    #process.hbhereco
    
#    process.particleFlowCluster
    process.ecalClusters*
    
    process.electronSeeds    #produced merged collection of TkDriven and Ecaldriven seeds
   # process.electronCkfTrackCandidates*process.electronGsfTracks #run electron tracking
    )
###############
outdir = "out_tests/"
outfilename = outdir + "reGsfTracking_maxCand_" + str(maxCand) + "_MaxChi2_" + str(maxChi2) + "_nSigma_" + str(nSigma) + ".root"
print "Writing output to file: " + outfilename

process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECODEBUGEventContent.outputCommands,
    fileName = cms.untracked.string( outfilename ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECODEBUG')
    )
)

process.RECODEBUGoutput.outputCommands.extend(cms.untracked.vstring('keep *_elGsfTracksWithQuality_*_*'))

# paths
process.p = cms.Path(
    process.myGsfReco 
)

process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)

process.schedule = cms.Schedule(
      process.p,
 #     process.RECODEBUGoutput_step
)


