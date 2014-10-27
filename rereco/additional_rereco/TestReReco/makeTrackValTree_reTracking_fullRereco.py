import FWCore.ParameterSet.Config as cms

process = cms.Process("reGsfTracking")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger = cms.Service("MessageLogger", #??
#                                    default = cms.untracked.PSet( limit = cms.untracked.int32(300) )
#                                    )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

source = cms.Source ("PoolSource",
                     fileNames=cms.untracked.vstring(
        #    'file:007CEDE1-B1D1-E311-9EC9-02163E00E9CC.root' 
        #    'file:step2.root'
        #    'file:SingleElectronPt10_RECO.root',
        #'/store/user/liis/GSF_tracking_samples/RelValZll_HLTDEBUG_PU50/007FA6E4-EC13-E411-A7A3-002618943866.root' #no evts(?)
        '/store/user/liis/GSF_tracking_samples/RelValZll_HLTDEBUG_PU50/00E5A57A-3D13-E411-99D5-0025905A6076.root'
        
#        'file:test2_sam.root' ## the last one
        #    'file:../EvtGeneration/SingleElectronPt10_RECO.root'
        #    'file:rawToReco.root'
        # -------- Zee produced by sam ----------------
        #        '/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_101_1_Fsb.root',
        #        '/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_100_2_qzC.root',
        #        '/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_102_1_VzM.root',
        #-----------------------------------------------
#        'file:SingleElectronPt10_RECO.root'   
        ),
#                     secondaryFileNames=secFiles #Provide corresponding DIGI files for tracking particles
                     )

process.source = source

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '') # dont use this one to rerun the reconstruction
#---------------------------------------
#An exception of category 'StdException' occurred while
#Exception Message:
#A std::exception was thrown.
#bitset::set
#---------------------------------------

print process.GlobalTag.globaltag #update automatically
#process.GlobalTag.globaltag = 'PRE_STA71_V4::All'
#process.GlobalTag.globaltag = 'START71_V8::All'


#process.load("Configuration.Geometry.GeometryIdeal_cff") #REMOVED
#process.load('Configuration.StandardSequences.GeometryPilot2_cff') #REMOVED
#process.load("Configuration.StandardSequences.MagneticField_cff") #remove

process.load("CondCore.DBCommon.CondDBSetup_cfi")

from SimGeneral.MixingModule.trackingTruthProducer_cfi import *
process.TrajectoryBuilderForElectrons.estimator = cms.string('Chi2') #comment out for an alternative trajectory finder
# 'Chi2A' -- separate costum producer defined in /TrackingTools/KalmanUpdators/python/Chi2MeasurementEstimatorESProducer_cfi.py
# TrajectoryBuilderForElectrons -- defined at TrackingTools/GsfTracking/python/CkfElectronCandidateMaker_cff.py: TrajectoryBuilderForElectrons =RecoTracker.CkfPattern.CkfTrajectoryBuilder_cfi.CkfTrajectoryBuilder.clone()

maxCandDefault = 5
maxChi2Default = 2000
nSigmaDefault = 3.0

maxCand = maxCandDefault
maxChi2 = 1 #maxChi2Default
nSigma = nSigmaDefault

########################################################################
# to change parameters  as in slides from A.Tropiano

process.TrajectoryBuilderForElectrons.maxCand = cms.int32( maxCand )
process.ElectronChi2.MaxChi2 = cms.double( maxChi2 )
process.ElectronChi2.nSigma = cms.double( nSigma )

########################################################################

process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi")

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.load("RecoTracker.FinalTrackSelectors.selectHighPurity_cfi")
process.elTracksWithQuality = process.selectHighPurity.clone(
    max_d0 = 0.02,
    minNumberLayers = 10,
#    src = "electronGsfTracks",
#    src = "electronGsfTracks",
)

#process.elGsfTracksWithQuality.src = cms.InputTag("electronGsfTracks")

# sequence for re-running gsfTracking over RECO
process.myGsfReco = cms.Sequence(
    #---- try to run on RelVal HLTDEBUG ----
  
    process.siPixelDigis
    *process.siStripDigis
    *process.ecalPreshowerDigis
    *process.ecalDigis
   
    *process.ecalGlobalUncalibRecHit
    *process.ecalDetIdToBeRecovered
    *process.ecalRecHitSequence
    *process.particleFlowRecHitECAL
    *process.particleFlowRecHitPS
    *process.ecalClustersNoPFBox

    # --- to run full localreco ---
    *process.hcalDigis
    *process.castorDigis
    *process.muonDTDigis ## maybe this is not needed
    *process.muonCSCDigis
    *process.muonRPCDigis
    *process.csctfDigis

    # ------------------------------

##    *process.offlineBeamSpot
##    *process.siStripZeroSuppression
##    *process.siStripClusters
##    *process.siPixelClusters  
 
##    *process.particleFlowClusterPS
##    *process.pfClusteringECAL #used to segfault due to wrong geometry


##########################
        
    #------------ standard ------------------
##    *process.siPixelRecHits
##   *process.siStripMatchedRecHits #make local hits
##    *process.MeasurementTrackerEvent #ADD
##    *process.siPixelClusterShapeCache # needed to add when moving from CMSSW_7_1_0_pre5 to pre7

    # --------------------
##    *process.PixelLayerTriplets
##    *process.pixelTracks
##    *process.pixelVertices
    #---------- stuff from iterTracking -----
##    *process.InitialStep
##    *process.DetachedTripletStep
##    *process.LowPtTripletStep
##    *process.PixelPairStep
##    *process.MixedTripletStep
##    *process.PixelLessStep
##    *process.TobTecStep
    #-------------------
#    *process.earlyGeneralTracks
#    *process.standalonemuontracking
#    *process.muonSeededStep
#   *process.muonSeededTracksInOut
#    *process.muonSeededTracksOutIn
#    *process.preDuplicateMergingGeneralTracks
 #   *process.duplicateTrackCandidates
 #   *process.mergedDuplicateTracks
    #------------------------
#    *process.generalTracks
  
    *process.reconstruction
#    *process.localreco
#    *process.muonlocalreco
#    *process.dtlocalreco
#    *process.standalonemuontracking
 #   *process.ancientMuonSeed
 #   *process.standAloneMuons
 #   *process.localreco
#    *process.iterTracking
#    *process.localreco

#    *process.electronGsfTracking

 #   *process.pixelPairStepClusters
#    *process.mixedTripletStep #Clusters
 #   *process.pixelLessStepClusters

    *process.electronSeedsSeq #ADD
    *process.electronGsfTracking
    *process.electronSeeds    #produced merged collection of TkDriven and Ecaldriven seeds
    *process.electronCkfTrackCandidates
    *process.electronGsfTracks #run electron tracking

#--------- Full electron reconstruction --------    

#    *process.particleFlowRecHitHCAL #  No need to do these are present, just need to get them in a single vector https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h

#    *process.pfTrack
#    *process.printEventContent #misses
#----    RefCore: A request to resolve a reference to a product of type 'std::vector<reco::PFRecHit>' with ProductID '2:921'
#----    can not be satisfied because the product cannot be found.
#----    Probably the branch containing the product is not stored in the input file.
    
#    *process.pfTrackElec
#    *process.ecalDrivenGsfElectronCores

#    *process.electronSequence
#    *process.ecalDrivenGsfElectronCores
#    *process.ecalDrivenGsfElectrons
#    *process.gsfElectrons
#-------------------------------------------------
 
#    *process.elTracksWithQuality
    
)

outdir = "out_tests/"
#outfilename = outdir + "reGsfTracking_maxCand_" + str(maxCand) + "_MaxChi2_" + str(maxChi2) + "_nSigma_" + str(nSigma) + ".root"
outfilename = "trackValTree_reTrk_1.root"
print "Writing output to file: " + outfilename

process.TFileService = cms.Service("TFileService", # if save
                                   fileName = cms.string(outfilename)
                                   )

#--- irrelevant here, defien both reco and sim in MakeTree/MakeTrackValTree/python/maketrackvaltree_cfi.py
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.TrackAssociatorByHits.SimToRecoDenominator = cms.string("reco") #Quality_SimToReco = shared hits/#reco(or #sim)
#process.TrackAssociatorByHits.Quality_SimToReco = cms.double(0.75)
#process.TrackAssociatorByHits.AbsoluteNumberOfHits = cms.bool(True)
#---------------- high purity selection of reco::Tracks---------------
process.load("PhysicsTools.RecoAlgos.recoTrackSelector_cfi")
process.cutsRecoTracksHp = process.recoTrackSelector.clone()
process.cutsRecoTracksHp.quality = cms.vstring("highPurity")
process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
process.cutsRecoTracksHp.maxAbsEta = cms.double(2.5)


process.ValidationSelectors = cms.Sequence(
    process.cutsRecoTracksHp
    )

#-----------------------------Filter Zee decays------------------
process.zeeFilter = cms.EDFilter("XtoFFbarFilter",
                        src = cms.InputTag("genParticles"),
                        idMotherX = cms.vint32(23),
                        idDaughterF = cms.vint32(11),
                        idMotherY = cms.vint32(),
                        idDaughterG = cms.vint32(),
                        )


#--------------------------- tree maker --------------------------
process.load("MakeTree.MakeTrackValTree.maketrackvaltree_cfi") # for writing output to a flat tree
process.trackValTreeMaker.isGSF = cms.bool(True)
process.trackValTreeMaker.leadingVertexOnly = cms.bool(False)

if process.trackValTreeMaker.isGSF:
    print "Running analysis on electron GSF tracks"
else:
    print "Running analysis on generalTracks"

if process.trackValTreeMaker.leadingVertexOnly:
    print "Require reco tracks to originate from the leading vertex"
else:
    print "Skip matching to leading vertex"

process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.preValidation = cms.Sequence(
    process.simHitTPAssocProducer
    )

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

# paths
process.p = cms.Path(
    process.zeeFilter
    *process.myGsfReco 
    *process.ValidationSelectors
#   ## *process.elTracksWithQuality #preselection for standard reco tracks
    *process.preValidation

#    *process.printEventContent 
    *process.trackValTreeMaker
    )


process.schedule = cms.Schedule(
    process.p
    )
