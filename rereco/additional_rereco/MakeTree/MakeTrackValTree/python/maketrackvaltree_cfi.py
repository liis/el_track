import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * 

TrackingParticleSelectionForEfficiency.tipTP = cms.double(3)
TrackingParticleSelectionForEfficiency.pdgIdTP = cms.vint32([-11, 11])

trackValTreeMaker = cms.EDAnalyzer('MakeTrackValTree',                                   
                                   TrackingParticleSelectionForEfficiency, # default tracking particle selection tresholds
                                   
                                   isGSF = cms.bool(True),
                                   trackLabelGSF = cms.InputTag("electronGsfTracks"),
                                   trackLabel = cms.InputTag("generalTracks"),
                                   elSeedLabel = cms.InputTag("electronMergedSeeds"),                      
)
