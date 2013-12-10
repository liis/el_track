import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * # default sim track selection tresholds

TrackingParticleSelectionForEfficiency.tipTP = cms.double(3)
TrackingParticleSelectionForEfficiency.pdgIdTP = cms.vint32([-11, 11])

TrackValTreeMaker = cms.EDAnalyzer('MakeTrackValTree',
                                   TrackingParticleSelectionForEfficiency,
                                   
                                   isGSF = cms.bool(True),
                                   recoTrackLabel = cms.InputTag()

                                   debug = cms.bool(False)
)
