import FWCore.ParameterSet.Config as cms
from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import * # default sim track selection tresholds
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import *

#from SimTracker.TrackAssociation.trackingParticleRecoTrackAssociation_cfi import * #trackingParticleRecoTrackAssociation
#from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import *

TrackingParticleSelectionForEfficiency.tipTP = cms.double(3)
TrackingParticleSelectionForEfficiency.pdgIdTP = cms.vint32([-11, 11])

#trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("electronGsfTracks")

TrackValTreeMaker = cms.EDAnalyzer('MakeTrackValTree',
                                   TrackingParticleSelectionForEfficiency,
                                   
                                   isGSF = cms.bool(True),
                                   debug = cms.bool(False)
)
