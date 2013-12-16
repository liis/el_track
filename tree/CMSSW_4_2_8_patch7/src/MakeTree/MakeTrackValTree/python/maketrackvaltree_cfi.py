import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * # default sim track selection tresholds
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import *

#from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
#from SimTracker.TrackAssociation.trackingParticleRecoTrackAssociation_cfi import * #trackingParticleRecoTrackAssociation
#from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import *

TrackingParticleSelectionForEfficiency.tipTP = cms.double(3)
TrackingParticleSelectionForEfficiency.pdgIdTP = cms.vint32([-11, 11])

#trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("electronGsfTracks")

TrackValTreeMaker = cms.EDAnalyzer('MakeTrackValTree',
                                   TrackingParticleSelectionForEfficiency,

                                   isGSF = cms.bool(True),
                                   debug = cms.bool(False),
                                   hitdebug = cms.bool(False),
                                   trackLabelGSF = cms.InputTag("electronGsfTracks"),
                                   trackLabel = cms.InputTag("generalTracks"),

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
                                   
                                   #associatePixel = cms.bool(False),
                                   #ROUList = cms.vstring('TrackerHitsTIBLowTof',
                                   #                      'TrackerHitsTIBHighTof',
                                   #                      'TrackerHitsTIDLowTof',
                                   #                      'TrackerHitsTIDHighTof',
                                   #                      'TrackerHitsTOBLowTof',
                                   #                      'TrackerHitsTOBHighTof',
                                   #                      'TrackerHitsTECLowTof',
                                   #                      'TrackerHitsTECHighTof',
                                   #                      'TrackerHitsPixelBarrelLowTof',
                                    #                     'TrackerHitsPixelBarrelHighTof',
                                    #                     'TrackerHitsPixelEndcapLowTof',
                                    #                     'TrackerHitsPixelEndcapHighTof'),
                                   #trajectoryInput = cms.string('generalTracks'),
                                   #associateRecoTracks = cms.bool(False),
                                   #   string trajectoryInput = "rsWithMaterialTracks"
                                   #associateStrip = cms.bool(True)
                                   
                                   )
