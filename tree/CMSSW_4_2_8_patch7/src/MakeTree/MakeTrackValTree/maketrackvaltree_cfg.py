import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackValTreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")

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

### validation-specific includes
process.load("DQMServices.Components.EDMtoMEConverter_cff")
#process.load("Validation.RecoTrack.cuts_cff")
#process.load("Validation.RecoTrack.MultiTrackValidator_cff")
#process.load("Validation.Configuration.postValidation_cff")
#process.load("Validation.RecoTrack.TrackValidation_cff") # defines track associator

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:aod.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('trackValTree.root')
                                   )


process.load("MakeTree.MakeTrackValTree.maketrackvaltree_cfi")

#--------------------Configure Track association-------------------------
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("electronGsfTracks")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string("reco") #Quality_SimToReco = shared hits/#reco(or #sim)
process.TrackAssociatorByHits.Quality_SimToReco = cms.double(0.75)
#process.TrackAssociatorByHits.AbsoluteNumberOfHits = cms.bool(True)

#---------------- high purity selection of reco::Tracks---------------
process.load("PhysicsTools.RecoAlgos.recoTrackSelector_cfi")
process.cutsRecoTracksHp = process.recoTrackSelector.clone()
process.cutsRecoTracksHp.quality = cms.vstring("highPurity")
process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
process.cutsRecoTracksHp.maxAbsEta = cms.double(2.5)
#process.cutsRecoTracksHp.src = cms.InputTag("electronGsfTracks") #cant apply for GSF tracks due to the assumtion of RecoTrackCollection (not View<reco::Track>)

process.ValidationSelectors = cms.Sequence(
    process.cutsRecoTracksHp
    )

process.AssociationMapProducers = cms.Sequence(
    process.trackingParticleRecoTrackAsssociation
    )

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.ValidationSelectors *
    process.AssociationMapProducers *

#    process.printEventContent *    # dump of event content after PAT-tuple production
    process.TrackValTreeMaker
    )

process.schedule = cms.Schedule(
    process.p
    )
