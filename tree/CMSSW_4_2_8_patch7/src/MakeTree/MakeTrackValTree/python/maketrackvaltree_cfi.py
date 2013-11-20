import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * # default sim track selection tresholds


demo = cms.EDAnalyzer('MakeTrackValTree',
                      TrackingParticleSelectionForEfficiency,
                      debug = cms.bool(True)
)
