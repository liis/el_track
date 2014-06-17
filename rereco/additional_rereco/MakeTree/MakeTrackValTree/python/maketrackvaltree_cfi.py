import FWCore.ParameterSet.Config as cms

trackValTreeMaker = cms.EDAnalyzer('MakeTrackValTree',
                                   isGSF = cms.bool(True),
                                   trackLabelGSF = cms.InputTag("electronGsfTracks"),
                                   trackLabel = cms.InputTag("generalTracks"),
                                   elSeedLabel = cms.InputTag("electronMergedSeeds"),
                      
)
