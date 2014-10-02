import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * 
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import *


TrackingParticleSelectionForEfficiency.lipTP = cms.double(300)
TrackingParticleSelectionForEfficiency.tipTP = cms.double(3) # omit secondary decays
TrackingParticleSelectionForEfficiency.signalOnlyTP = cms.bool(True) # if True, dont consider PU tracks

TrackingParticleSelectionForEfficiency.minHitTP = cms.int32(3)
TrackingParticleSelectionForEfficiency.ptMinTP = cms.double(0.9)
TrackingParticleSelectionForEfficiency.pdgIdTP = cms.vint32([-11, 11])


import SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi 
TrackAssociatorByHitsRecoDenom= SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone(
    ComponentName = cms.string('TrackAssociatorByHitsRecoDenom'),  
    )

trackValTreeMaker = cms.EDAnalyzer('MakeTrackValTree',                                   
                                   TrackingParticleSelectionForEfficiency, # default tracking particle selection tresholds
                                   
                                   isGSF = cms.bool(True),
                                   isSinglePart = cms.bool(False),
                                   trackLabelGSF = cms.InputTag("electronGsfTracks"),
                                   trackLabel = cms.InputTag("generalTracks"),
                                   elSeedLabel = cms.InputTag("electronMergedSeeds"),                      

                                   #------------ for hitAssociation -----------------
                                   Quality_SimToReco = cms.double(0.5),
                                   associateRecoTracks = cms.bool(True),
                                   UseGrouped = cms.bool(True),
                                   associatePixel = cms.bool(True),
                                   ROUList = cms.vstring('TrackerHitsTIBLowTof',
                                                         'TrackerHitsTIBHighTof',
                                                         'TrackerHitsTIDLowTof',
                                                         'TrackerHitsTIDHighTof',
                                                         'TrackerHitsTOBLowTof',
                                                         'TrackerHitsTOBHighTof',
                                                         'TrackerHitsTECLowTof',
                                                         'TrackerHitsTECHighTof',
                                                         'TrackerHitsPixelBarrelLowTof',
                                                         'TrackerHitsPixelBarrelHighTof',
                                                         'TrackerHitsPixelEndcapLowTof',
                                                         'TrackerHitsPixelEndcapHighTof'),
                                   UseSplitting = cms.bool(True),
                                   ComponentName = cms.string('TrackAssociatorByHits'),
                                   UsePixels = cms.bool(True),
                                   ThreeHitTracksAreSpecial = cms.bool(True),
                                   AbsoluteNumberOfHits = cms.bool(False),
                                   associateStrip = cms.bool(True),
                                   Purity_SimToReco = cms.double(0.75),
                                   Cut_RecoToSim = cms.double(0.75),
                                   SimToRecoDenominator = cms.string('sim') ##"reco"
                                   
                                   )
