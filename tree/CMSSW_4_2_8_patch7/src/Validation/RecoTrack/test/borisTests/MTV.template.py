import FWCore.ParameterSet.Config as cms

process = cms.Process("MULTITRACKVALIDATOR")

# message logger
process.MessageLogger = cms.Service("MessageLogger",
     default = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

# source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( ['file:./step2.root'])

process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("Validation.RecoTrack.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("Validation.Configuration.postValidation_cff")
process.load("Validation.RecoTrack.TrackValidation_cff")

process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
#--- To change the fraction of sHits that are required to be matched by the associator
#    The default is 0.75
#process.TrackAssociatorByHits.Purity_SimToReco = cms.double(0.60)
#process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.60)
#---

########### configuration MultiTrackValidator ########
process.multiTrackValidator.outputFile = 'multitrackvalidator.root'
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.skipHistoFit=cms.untracked.bool(False)
process.multiTrackValidator.runStandalone=cms.bool(True)
#process.multiTrackValidator.label=cms.VInputTag(cms.InputTag("generalTracks"),
process.multiTrackValidator.label=cms.VInputTag(
                                   cms.InputTag("cutsRecoTracksHp"),
#                                   cms.InputTag("cutsRecoTracksZero"),
#                                   cms.InputTag("cutsRecoTracksZeroHp"),
#                                   cms.InputTag("cutsRecoTracksFirst"),
#                                   cms.InputTag("cutsRecoTracksFirstHp"),
#                                   cms.InputTag("cutsRecoTracksSecond"),
#                                   cms.InputTag("cutsRecoTracksSecondHp"),
#                                   cms.InputTag("cutsRecoTracksThird"),
#                                   cms.InputTag("cutsRecoTracksThirdHp"),
#                                   cms.InputTag("cutsRecoTracksFourth"),
#                                   cms.InputTag("cutsRecoTracksFourthHp"),
#                                   cms.InputTag("cutsRecoTracksFifth"),
#                                   cms.InputTag("cutsRecoTracksFifthHp")
                                   )
process.multiTrackValidator.useLogPt=cms.untracked.bool(True)
process.multiTrackValidator.minPt = cms.double(0.1)
process.multiTrackValidator.maxPt = cms.double(300.0)
process.multiTrackValidator.nintPt = cms.int32(40)

#--- Filter on TrackingParticles
#    pt in [0,2.8] when calculating the tracking Fake rate
#    pt in [0,2.5] when calculating the efficiency vs eta
#    pt in eta slice when calculating the efficiency vs pt for barrel/transition/endcap
process.multiTrackValidator.useLogPt=cms.untracked.bool(True)
process.multiTrackValidator.histoProducerAlgoBlock.minPt = cms.double(0.1)
process.multiTrackValidator.histoProducerAlgoBlock.maxPt = cms.double(300.0)
process.multiTrackValidator.histoProducerAlgoBlock.nintPt = cms.int32(40)
#CFG_TpBarrEff#process.multiTrackValidator.minAbsEtaTP = cms.double(0.0)
#CFG_TpBarrEff#process.multiTrackValidator.maxAbsEtaTP = cms.double(0.9)
#CFG_TpTraEff#process.multiTrackValidator.minAbsEtaTP = cms.double(0.9)
#CFG_TpTraEff#process.multiTrackValidator.maxAbsEtaTP = cms.double(1.4)
#CFG_TpEndEff#process.multiTrackValidator.minAbsEtaTP = cms.double(1.4)
#CFG_TpEndEff#process.multiTrackValidator.maxAbsEtaTP = cms.double(2.5)
#CFG_TpFakes#process.multiTrackValidator.minAbsEtaTP = cms.double(0.0)
#CFG_TpFakes#process.multiTrackValidator.maxAbsEtaTP = cms.double(2.8)
#CFG_TpVsEta#process.multiTrackValidator.minAbsEtaTP = cms.double(0.0)
#CFG_TpVsEta#process.multiTrackValidator.maxAbsEtaTP = cms.double(2.8)
#process.multiTrackValidator.minAbsEtaTP = cms.double(0.0)
#process.multiTrackValidator.maxAbsEtaTP = cms.double(2.5)
#---


#--- uncomment this part to run MTV on GsfTrack collection
#
#CFG_GSF#process.cutsRecoTracksHp = cms.EDFilter("RecoGsfTrackSelector",
#CFG_GSF#    src = cms.InputTag("electronGsfTracks"),
###    src = cms.InputTag("elGsfTracksWithQuality"),
#CFG_GSF#    beamSpot = cms.InputTag("offlineBeamSpot"),
#CFG_GSF#    algorithm = cms.vstring(), 
#CFG_GSF#    maxChi2 = cms.double(10000.0),
###    #quality = cms.vstring('highPurity'), ## NEW
###    quality = cms.vstring('loose'), ## NEW
#CFG_GSF#    quality = cms.vstring(), ## NEW
#CFG_GSF#    minRapidity = cms.double(-5.0),
#CFG_GSF#    maxRapidity = cms.double(5.0),
#CFG_GSF#    tip = cms.double(120.0),
#CFG_GSF#    lip = cms.double(300.0),
#CFG_GSF#    ptMin = cms.double(0.1),
#CFG_GSF#    min3DHit = cms.int32(0),
#CFG_GSF#    minHit = cms.int32(3),
#CFG_GSF#    minAbsEta = cms.double(1.4),
#CFG_GSF#    maxAbsEta = cms.double(2.5)
#CFG_GSF#)
#CFG_GSF#process.multiTrackValidator.histoProducerAlgoBlock.useGsf = cms.bool(True)
#---

#--- Filter on track collection
#    pt in [0,2.8] when calculating the tracking efficiency
#    pt in eta slice when calculating the fake rate vs pt for barrel/transition/endcap
#CFG_TkBarrFake#process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
#CFG_TkBarrFake#process.cutsRecoTracksHp.maxAbsEta = cms.double(0.9)
#CFG_TkTraFake#process.cutsRecoTracksHp.minAbsEta = cms.double(0.9)
#CFG_TkTraFake#process.cutsRecoTracksHp.maxAbsEta = cms.double(1.4)
#CFG_TkEndFake#process.cutsRecoTracksHp.minAbsEta = cms.double(1.4)
#CFG_TkEndFake#process.cutsRecoTracksHp.maxAbsEta = cms.double(2.5)
#CFG_TkVsEta#process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
#CFG_TkVsEta#process.cutsRecoTracksHp.maxAbsEta = cms.double(2.8)
#CFG_TkEff#process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
#CFG_TkEff#process.cutsRecoTracksHp.maxAbsEta = cms.double(2.8)
#process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
#process.cutsRecoTracksHp.maxAbsEta = cms.double(2.5)



process.multiTrackValidator.UseAssociators = cms.bool(True)
process.ValidationSelectors = cms.Sequence( process.cutsRecoTracksHp
#                                process.cutsRecoTracksZero*
#                                process.cutsRecoTracksZeroHp*
#                                process.cutsRecoTracksFirst*
#                                process.cutsRecoTracksFirstHp*
#                                process.cutsRecoTracksSecond*
#                                process.cutsRecoTracksSecondHp*
#                                process.cutsRecoTracksThird*
#                                process.cutsRecoTracksThirdHp*
#                                process.cutsRecoTracksFourth*
#                                process.cutsRecoTracksFourthHp*
#                                process.cutsRecoTracksFifth*
#                                process.cutsRecoTracksFifthHp 
)

process.validation = cms.Sequence(
    process.multiTrackValidator
)



# paths
process.p = cms.Path(
    process.ValidationSelectors *
    process.validation
)
process.schedule = cms.Schedule(
      process.p
)


