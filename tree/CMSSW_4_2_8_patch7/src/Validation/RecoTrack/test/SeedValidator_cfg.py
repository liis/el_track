import FWCore.ParameterSet.Config as cms

process = cms.Process("SEEDVALIDATOR")
process.load("Configuration/StandardSequences/GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START42_V11::All' #'MC_31X_V2::All'
#process.MessageLogger.categories = ['TrackAssociator', 'TrackValidator']
#process.MessageLogger.debugModules = ['*']
#process.MessageLogger.cout = cms.untracked.PSet(
#    threshold = cms.untracked.string('DEBUG'),
#    default = cms.untracked.PSet(
#        limit = cms.untracked.int32(0)
#    ),
#    TrackAssociator = cms.untracked.PSet(
#        limit = cms.untracked.int32(0)
#    ),
#    TrackValidator = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#    )
#)
#process.MessageLogger.cerr = cms.untracked.PSet(
#    placeholder = cms.untracked.bool(True)
#)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring([
    'file:borisTests/aod.root'
]                                     ),
   secondaryFileNames=cms.untracked.vstring([
]
))

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("Validation.RecoTrack.cuts_cff")

process.load("Validation.RecoTrack.TrackerSeedValidator_cff")
#process.multiTrackValidator.associators = cms.vstring('TrackAssociatorByHits','TrackAssociatorByChi2')
#process.multiTrackValidator.UseAssociators = True
#process.multiTrackValidator.label = ['cutsRecoTracks']
#process.multiTrackValidator.label_tp_effic = cms.InputTag("cutsTPEffic")
#process.multiTrackValidator.label_tp_fake  = cms.InputTag("cutsTPFake")
#process.multiTrackValidator.associatormap = cms.InputTag(assoc2GsfTracks)
process.trackerSeedValidator.outputFile = 'file.root'

# Tracking Truth and mixing module, if needed
#process.load("SimGeneral.MixingModule.mixNoPU_cfi")
#process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")

process.evtInfo = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits*process.ckftracks*process.cutsTPEffic*process.cutsTPFake*process.trackerSeedValidator)
process.p = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits*process.ckftracks*process.trackerSeedValidator)
#process.p = cms.Path(process.multiTrackValidator)
process.ep = cms.EndPath(process.evtInfo)


