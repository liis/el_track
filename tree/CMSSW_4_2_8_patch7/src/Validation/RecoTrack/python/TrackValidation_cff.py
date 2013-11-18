import FWCore.ParameterSet.Config as cms

import SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi 
import SimTracker.TrackAssociation.TrackAssociatorByHits_cfi 
import Validation.RecoTrack.MultiTrackValidator_cfi
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import *
from SimTracker.TrackAssociation.CosmicParametersDefinerForTP_cfi import *
from Validation.RecoTrack.PostProcessorTracker_cfi import *
import PhysicsTools.RecoAlgos.recoTrackSelector_cfi

TrackAssociatorByHitsRecoDenom= SimTracker.TrackAssociation.TrackAssociatorByHits_cfi.TrackAssociatorByHits.clone(
    ComponentName = cms.string('TrackAssociatorByHitsRecoDenom'),  
    SimToRecoDenominator = cms.string('reco')
    )
# Validation iterative steps
cutsRecoTracksZero = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksZero.algorithm=cms.vstring("iter0")

cutsRecoTracksFirst = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksFirst.algorithm=cms.vstring("iter1")

cutsRecoTracksSecond = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksSecond.algorithm=cms.vstring("iter2")

cutsRecoTracksThird = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksThird.algorithm=cms.vstring("iter3")

cutsRecoTracksFourth = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksFourth.algorithm=cms.vstring("iter4")

cutsRecoTracksFifth = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksFifth.algorithm=cms.vstring("iter5")

# high purity
cutsRecoTracksHp = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksHp.quality=cms.vstring("highPurity")

cutsRecoTracksZeroHp = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksZeroHp.algorithm=cms.vstring("iter0")
cutsRecoTracksZeroHp.quality=cms.vstring("highPurity")

cutsRecoTracksFirstHp = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksFirstHp.algorithm=cms.vstring("iter1")
cutsRecoTracksFirstHp.quality=cms.vstring("highPurity")

cutsRecoTracksSecondHp = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksSecondHp.algorithm=cms.vstring("iter2")
cutsRecoTracksSecondHp.quality=cms.vstring("highPurity")

cutsRecoTracksThirdHp = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksThirdHp.algorithm=cms.vstring("iter3")
cutsRecoTracksThirdHp.quality=cms.vstring("highPurity")

cutsRecoTracksFourthHp = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksFourthHp.algorithm=cms.vstring("iter4")
cutsRecoTracksFourthHp.quality=cms.vstring("highPurity")

cutsRecoTracksFifthHp = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
cutsRecoTracksFifthHp.algorithm=cms.vstring("iter5")
cutsRecoTracksFifthHp.quality=cms.vstring("highPurity")

trackValidator= Validation.RecoTrack.MultiTrackValidator_cfi.multiTrackValidator.clone()

trackValidator.label=cms.VInputTag(cms.InputTag("generalTracks"),
                                   cms.InputTag("cutsRecoTracksHp"),
                                   cms.InputTag("cutsRecoTracksZero"),
                                   cms.InputTag("cutsRecoTracksZeroHp"),
                                   cms.InputTag("cutsRecoTracksFirst"),
                                   cms.InputTag("cutsRecoTracksFirstHp"),
                                   cms.InputTag("cutsRecoTracksSecond"),
                                   cms.InputTag("cutsRecoTracksSecondHp"),
                                   cms.InputTag("cutsRecoTracksThird"),
                                   cms.InputTag("cutsRecoTracksThirdHp"),
                                   cms.InputTag("cutsRecoTracksFourth"),
                                   cms.InputTag("cutsRecoTracksFourthHp"),
                                   cms.InputTag("cutsRecoTracksFifth"),
                                   cms.InputTag("cutsRecoTracksFifthHp")
                                   )
trackValidator.skipHistoFit=cms.untracked.bool(True)
trackValidator.useLogPt=cms.untracked.bool(True)
#trackValidator.minpT = cms.double(-1)
#trackValidator.maxpT = cms.double(3)
#trackValidator.nintpT = cms.int32(40)

# the track selectors
tracksValidationSelectors = cms.Sequence( cutsRecoTracksHp*
                                cutsRecoTracksZero*
                                cutsRecoTracksZeroHp*
                                cutsRecoTracksFirst*
                                cutsRecoTracksFirstHp*
                                cutsRecoTracksSecond*
                                cutsRecoTracksSecondHp*
                                cutsRecoTracksThird*
                                cutsRecoTracksThirdHp*
                                cutsRecoTracksFourth*
                                cutsRecoTracksFourthHp*
                                cutsRecoTracksFifth*
                                cutsRecoTracksFifthHp )

# selectors go into separate "prevalidation" sequence
tracksValidation = cms.Sequence( trackValidator)
tracksValidationFS = cms.Sequence( trackValidator )

