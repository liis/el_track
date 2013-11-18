import FWCore.ParameterSet.Config as cms

process = cms.Process("reTracking")

# message logger
process.MessageLogger = cms.Service("MessageLogger",
     default = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

# source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
#readFiles.extend( ['root://pcmssd12//data/mangano/SingleElMinusPt10.root'])
readFiles.extend( ['root://pcmssd12//data/mangano/SingleElFlatPt.root'])

process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

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




#### re-reco of gsfElectrons ####
process.electronCkfTrackCandidates.src = cms.InputTag("newCombinedSeeds")
process.newCombinedSeeds.seedCollections = cms.VInputTag(
     cms.InputTag('newSeedFromTriplets'),
     cms.InputTag('newSeedFromPairs'),
     cms.InputTag('secTriplets'),
     cms.InputTag('thTriplets'), 
     cms.InputTag('fourthPLSeeds'),
)


process.elGsfTracksQualityLoose = cms.EDProducer("AnalyticalGsfTrackSelector",
    max_d0 = cms.double(100.0),
    minNumber3DLayers = cms.uint32(0),
    applyAbsCutsIfNoPV = cms.bool(False),
    qualityBit = cms.string('loose'),
    minNumberLayers = cms.uint32(0),
    chi2n_par = cms.double(1.6),
    nSigmaZ = cms.double(3.0),
    dz_par2 = cms.vdouble(0.45, 4.0),
    applyAdaptedPVCuts = cms.bool(True),
    dz_par1 = cms.vdouble(0.65, 4.0),
    copyTrajectories = cms.untracked.bool(True),
    vtxNumber = cms.int32(-1),
    keepAllTracks = cms.bool(False),
    maxNumberLostLayers = cms.uint32(999),
    beamspot = cms.InputTag("offlineBeamSpot"),
    copyExtras = cms.untracked.bool(True),
    vertexCut = cms.string('ndof>=2&!isFake'),
    max_z0 = cms.double(100.0),
    src = cms.InputTag("preFilterZeroStepTracks"),
    vertices = cms.InputTag("pixelVertices"),
    d0_par2 = cms.vdouble(0.55, 4.0),
    d0_par1 = cms.vdouble(0.55, 4.0),
    res_par = cms.vdouble(0.003, 0.01)
)

process.elGsfTracksWithQuality = cms.EDProducer("AnalyticalGsfTrackSelector",
    max_d0 = cms.double(100.0),
    minNumber3DLayers = cms.uint32(3),
    applyAbsCutsIfNoPV = cms.bool(False),
    qualityBit = cms.string('highPurity'),
    minNumberLayers = cms.uint32(3),
    chi2n_par = cms.double(0.7),
    nSigmaZ = cms.double(3.0),
    dz_par2 = cms.vdouble(0.4, 4.0),
    applyAdaptedPVCuts = cms.bool(True),
    dz_par1 = cms.vdouble(0.35, 4.0),
    copyTrajectories = cms.untracked.bool(True),
    vtxNumber = cms.int32(-1),
    keepAllTracks = cms.bool(True),
    maxNumberLostLayers = cms.uint32(2),
    beamspot = cms.InputTag("offlineBeamSpot"),
    copyExtras = cms.untracked.bool(True),
    vertexCut = cms.string('ndof>=2&!isFake'),
    max_z0 = cms.double(100.0),
    src = cms.InputTag("zeroStepWithTightQuality"),
    vertices = cms.InputTag("pixelVertices"),
    d0_par2 = cms.vdouble(0.4, 4.0),
    d0_par1 = cms.vdouble(0.3, 4.0),
    res_par = cms.vdouble(0.003, 0.001)
)

process.elGsfTracksQualityLoose.src = cms.InputTag("electronGsfTracks")
process.elGsfTracksWithQuality.src = cms.InputTag("elGsfTracksQualityLoose")

process.newCombinedSeeds.clusterRemovalInfos = cms.VInputTag(
     cms.InputTag(''),
     process.preFilterStepOneTracks.clusterRemovalInfo,
     process.secWithMaterialTracks.clusterRemovalInfo,
     process.thWithMaterialTracks.clusterRemovalInfo,
     process.fourthWithMaterialTracks.clusterRemovalInfo
     )

process.myGsfReco = cms.Sequence(
    process.siPixelRecHits+process.siStripMatchedRecHits+
    process.iterTracking+
    process.newCombinedSeeds+
    process.electronCkfTrackCandidates+process.electronGsfTracks+
    process.elGsfTracksQualityLoose+
    process.elGsfTracksWithQuality)
###############


process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECODEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('step2.root'),
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
      process.p,process.RECODEBUGoutput_step
)


