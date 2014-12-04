import FWCore.ParameterSet.Config as cms

process = cms.Process("reGsfTracking")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

#process.Timing = cms.Service("Timing")

# source
readFiles = cms.untracked.vstring()
#------------------------- define secondary files -----------------------------
useSecFiles = False #needed for MakeTrackValTree (tracking particles)
run_PSI     = False # only relevant if using sec files

if not useSecFiles: #Use secondary files
    secFiles = cms.untracked.vstring()
elif run_PSI:
    secFiles = cms.untracked.vstring(
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/02A8EF51-62D1-E311-9270-02163E00EB1C.root',
        )
else:  # run at EE
    secFiles = scms.untracked.vstring(
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/02A8EF51-62D1-E311-9270-02163E00EB1C.root',
        )
#----------------------------------------------------------------------------------------------

#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
    
source = cms.Source ("PoolSource",
                     fileNames=cms.untracked.vstring(
        #    'file:007CEDE1-B1D1-E311-9EC9-02163E00E9CC.root' 
        #    'file:step2.root'
        #    'file:SingleElectronPt10_RECO.root',
        
        ########### work in progress ##############
        #    '/store/user/liis/GSF_tracking_samples/RelValZll_HLTDEBUG_PU50/007FA6E4-EC13-E411-A7A3-002618943866.root'
        #    '/store/user/liis/GSF_tracking_samples/RelValZll_HLTDEBUG_PU50/00E5A57A-3D13-E411-99D5-0025905A6076.root'
        ###########################################
        
#        'file:../EvtGeneration/file:test2_sam.root' ## the last one
        'file:../EvtGeneration/zee_reco_all.root'
#        'file:../EvtGeneration/zee_reco_standard.root'
#        'file:../EvtGeneration/zee_reco.root'
        
        #        'file:../EvtGeneration/samtest_reco.root'
        #    'file:rawToReco.root'
        # -------- Zee produced by sam ----------------


        #'/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_101_1_Fsb.root',
        #'/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_100_2_qzC.root',
        #'/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_102_1_VzM.root',
        #-----------------------------------------------
        #        'file:SingleElectronPt10_RECO.root'   
        ),
                     secondaryFileNames=secFiles #Provide corresponding DIGI files for tracking particles
                     )

process.source = source

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '') ## dont use
print "Global tag = " + str(process.GlobalTag.globaltag)

process.load("CondCore.DBCommon.CondDBSetup_cfi")

from SimGeneral.MixingModule.trackingTruthProducer_cfi import *

#process.TrajectoryBuilderForElectrons.estimator = cms.string('Chi2') #comment out for an alternative trajectory finder
# 'Chi2A' -- separate costum producer defined in /TrackingTools/KalmanUpdators/python/Chi2MeasurementEstimatorESProducer_cfi.py
# TrajectoryBuilderForElectrons -- defined at TrackingTools/GsfTracking/python/CkfElectronCandidateMaker_cff.py: TrajectoryBuilderForElectrons =RecoTracker.CkfPattern.CkfTrajectoryBuilder_cfi.CkfTrajectoryBuilder.clone()
## ---- !!! if this is uncommented, normal overwriting of tracking parameters doesn't work !!! -----


maxCandDefault = 5
maxChi2Default = 2000
nSigmaDefault = 3.0

maxCand = maxCandDefault
maxChi2 = 20 #maxChi2Default
nSigma = nSigmaDefault

########################################################################
# to change parameters  as in slides from A.Tropiano

process.TrajectoryBuilderForElectrons.maxCand = cms.int32( maxCand )
process.ElectronChi2.MaxChi2 = cms.double( maxChi2 )
process.ElectronChi2.nSigma = cms.double( nSigma )

########################################################################

process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi") #?

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.load("RecoTracker.FinalTrackSelectors.selectHighPurity_cfi")
process.elTracksWithQuality = process.selectHighPurity.clone(
    max_d0 = 0.02,
    minNumberLayers = 10,
#    src = "electronGsfTracks",
#    src = "electronGsfTracks",
)

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.preGsfReco = cms.Sequence(
    process.siPixelDigis
    *process.siPixelClusters
    *process.siStripDigis
    *process.siStripZeroSuppression
    )

process.additionalReco = cms.Sequence(
    process.siPixelRecHits
    *process.siStripMatchedRecHits
    )

# sequence for re-running gsfTracking over RECO
process.myGsfReco = cms.Sequence(
    process.siPixelRecHits
    *process.siStripMatchedRecHits #make local hits
    *process.MeasurementTrackerEvent #ADD
    *process.siPixelClusterShapeCache # needed to add when moving from CMSSW_7_1_0_pre5 to pre7

    *process.iterTracking #probably don't need to run all of it
    
    *process.electronSeedsSeq #ADD
    *process.electronSeeds    #produced merged collection of TkDriven and Ecaldriven seeds
    *process.electronCkfTrackCandidates
    *process.electronGsfTracks #run electron tracking

    # --------- for full electron reconstruction ----------
    # --> note that this doesn't work on sam Zee sample due to some missing collections of type 'std::vector<reco::PFRecHit>' 

#    *process.pfTrack #contained in pfTrackingGlobalReco
#    *process.pfTrackElec
    #--------- for particleFlowTmp producing pfCandidates ----------
    
    *process.globalMuons
    *process.muons1stStep
##    *process.particleFlowCluster # <-- upper 3 for running particleFlowBlock
    

    *process.pfTrackingGlobalReco

    *process.ecalDrivenGsfElectronCores #after pfTrackingGlobalReco, before electronsWithPresel
    *process.ecalDrivenGsfElectrons

##    *process.electronsWithPresel
##    *process.mvaElectrons


#    *process.particleFlowBlock #attempt to get PFCandidateEGammaExtra for, but not helping
#    *process.particleFlowEGamma

#    *process.gedPhotonsTmp
#    *process.gedGsfElectronCores
#    *process.gedGsfElectronsTmp
 
##    *process.particleFlowReco ## rereco everything -- works if full output is written out
#    *process.particleFlowClusterECALUncorrected
#    *process.particleFlowClusterECAL
#    *process.particleFlowEGammaFull
    *process.particleFlowTmp
    *process.particleFlowTmpPtrs
#    *process.particleFlowEGammaFinal

#
    *process.pfSelectedElectrons
#    *process.pfElectronTranslator #Sequence # maybe full sequence is not needed

#    *process.gsfElectronCores
#    *process.gsfElectrons
#    *process.gsfElectronSequence
    #-------------------------
    )

outdir = "out_tests/"
outfilename = "trackValTree_reTrk_newtry2.root"
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

#--------------------------- tree maker --------------------------
process.load("MakeTree.MakeTrackValTree.maketrackvaltree_cfi") # for writing output to a flat tree
process.trackValTreeMaker.isGSF = cms.bool(True)
process.trackValTreeMaker.leadingVertexOnly = cms.bool(False) # consider only tracks from the leading vertex (needed for Zee sample without the PU tracking particles)

#-----------------------------Filter Zee decays------------------                                                                       
process.zeeFilter = cms.EDFilter("XtoFFbarFilter",
                                 src = cms.InputTag("genParticles"),
                                 idMotherX = cms.vint32(23),
                                 idDaughterF = cms.vint32(11),
                                 idMotherY = cms.vint32(),
                                 idDaughterG = cms.vint32(),
                                 )
#---------------------------------------------------------------


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

# paths
process.p = cms.Path(
    process.zeeFilter # Sams 13 TeV sample is already filtered for electron decays
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
