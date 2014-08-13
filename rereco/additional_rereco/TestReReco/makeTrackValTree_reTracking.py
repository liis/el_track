import FWCore.ParameterSet.Config as cms

process = cms.Process("reGsfTracking")

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger = cms.Service("MessageLogger", #??
                                    default = cms.untracked.PSet( limit = cms.untracked.int32(100) )
                                    )

# source
readFiles = cms.untracked.vstring()
#------------------------- define secondary files -----------------------------
useSecFiles = False #needed for MakeTrackValTree (tracking particles)
run_PSI     = True

if not useSecFiles: #Use secondary files
    secFiles = cms.untracked.vstring()
elif run_PSI:
    secFiles = cms.untracked.vstring(
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/02A8EF51-62D1-E311-9270-02163E00EB1C.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/22360CB9-62D1-E311-AAC8-02163E00E8F4.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/2C4CD1E2-62D1-E311-840A-02163E00EA17.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/2CC8BC3A-62D1-E311-80C4-02163E00EAC3.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/3E92DFF4-62D1-E311-B169-02163E00E83C.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/4AF6D3E5-62D1-E311-BB13-02163E00F335.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/4E9486E1-62D1-E311-B091-02163E00EB07.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/563694CA-62D1-E311-B8B7-02163E00E9F6.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/6A17CF7C-62D1-E311-A943-02163E00E6D6.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/76B6063A-72D1-E311-95E6-02163E00CDE3.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/9E65F9A2-62D1-E311-8D37-02163E00E8D2.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/9EA489AB-62D1-E311-9A1A-02163E00E776.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/AE3BAC91-62D1-E311-B275-02163E00E869.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/C637DD05-63D1-E311-A4FB-02163E00EAC7.root',
        '/store/user/liis/GSF_tracking_samples/Zee_DIGI_PU50_710pre7/DE0B8342-62D1-E311-AECF-02163E00EA85.root'
        )
else:  # run at EE
    secFiles = scms.untracked.vstring(
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/02A8EF51-62D1-E311-9270-02163E00EB1C.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/22360CB9-62D1-E311-AAC8-02163E00E8F4.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/2C4CD1E2-62D1-E311-840A-02163E00EA17.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/2CC8BC3A-62D1-E311-80C4-02163E00EAC3.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/3E92DFF4-62D1-E311-B169-02163E00E83C.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/4AF6D3E5-62D1-E311-BB13-02163E00F335.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/4E9486E1-62D1-E311-B091-02163E00EB07.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/563694CA-62D1-E311-B8B7-02163E00E9F6.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/6A17CF7C-62D1-E311-A943-02163E00E6D6.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/76B6063A-72D1-E311-95E6-02163E00CDE3.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/9E65F9A2-62D1-E311-8D37-02163E00E8D2.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/9EA489AB-62D1-E311-9A1A-02163E00E776.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/AE3BAC91-62D1-E311-B275-02163E00E869.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/C637DD05-63D1-E311-A4FB-02163E00EAC7.root',
        'file:/hdfs/cms/store/user/liis/El_GSF_studies/Zee_7_1_0_pre7_PU50_DIGI/DE0B8342-62D1-E311-AECF-02163E00EA85.root', 
        )
#----------------------------------------------------------------------------------------------
    
source = cms.Source ("PoolSource",
                     fileNames=cms.untracked.vstring(
#    'file:007CEDE1-B1D1-E311-9EC9-02163E00E9CC.root' 
#    'file:step2.root'
#    'file:samtest_reco.root'
#    'file:test2_sam.root' ## the last one
#    'file:../EvtGeneration/outfiles/testing/SingleElectronPt10_RECO.root'
#    'file:test_sam_zee.root'
#    'file:rawToReco.root'
# -------- Zee produced by sam ----------------
#        '/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_101_1_Fsb.root',
#        '/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_100_2_qzC.root',
#        '/store/group/phys_egamma/sharper/DYJetsToLL_M-50_13TeV-pythia6/EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1/ffac44eb0cb582bdcc6ecfb3c5f327a8/DYJetsToLL_M-50_13TeV-pythia6_EGM711_PU40bx25_POSTLS171_V11-v1_102_1_VzM.root',
#-----------------------------------------------
        'file:SingleElectronPt10_RECO.root'   
        ),
                     secondaryFileNames=secFiles #Provide corresponding DIGI files for tracking particles
                     )

process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
process.GlobalTag.globaltag = 'PRE_STA71_V4::All'

process.load('Configuration/StandardSequences/Services_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")

from SimGeneral.MixingModule.trackingTruthProducer_cfi import *

#process.TrajectoryBuilderForElectrons.estimator = cms.string('Chi4A')

maxCandDefault = 5
maxChi2Default = 2000
nSigmaDefault = 3.0

maxCand = 6
maxChi2 = 100
nSigma = 4

########################################################################
# to change parameters  as in slides from A.Tropiano

process.TrajectoryBuilderForElectrons.maxCand = cms.int32( maxCand )
process.ElectronChi2.MaxChi2 = cms.double( maxChi2 )
process.ElectronChi2.nSigma = cms.double( nSigma )

########################################################################

process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi")

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.elGsfTracksWithQuality = cms.EDProducer("AnalyticalTrackSelector",
    src = cms.InputTag("generalTracks"),
    keepAllTracks = cms.bool(False), ## if set to true tracks failing this filter are kept in the output
    beamspot = cms.InputTag("offlineBeamSpot"),

    # vertex selection 
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("pixelVertices"),
    vtxNumber = cms.int32(-1),
    vertexCut = cms.string('ndof>=2&!isFake'),

    #untracked bool copyTrajectories = true // when doing retracking before
    copyTrajectories = cms.untracked.bool(False),
    copyExtras = cms.untracked.bool(True), ## set to false on AOD
    qualityBit = cms.string('tight'), ## set to '' or comment out if you don't want to set the bit

    # parameters for adapted optimal cuts on chi2 and primary vertex compatibility
    chi2n_par = cms.double(0.7),
    chi2n_no1Dmod_par = cms.double(9999.),
    res_par = cms.vdouble(0.003, 0.01),
    d0_par1 = cms.vdouble(0.3, 4.0),
    dz_par1 = cms.vdouble(0.35, 4.0),
    d0_par2 = cms.vdouble(0.4, 4.0),
    dz_par2 = cms.vdouble(0.4, 4.0),
    # Boolean indicating if adapted primary vertex compatibility cuts are to be applied.
    applyAdaptedPVCuts = cms.bool(True),

    # Impact parameter absolute cuts.
    max_d0 = cms.double(100.),
    max_z0 = cms.double(100.),
    nSigmaZ = cms.double(4.),

    # Cuts on numbers of layers with hits/3D hits/lost hits. 
    minNumberLayers = cms.uint32(3),
    minNumber3DLayers = cms.uint32(3),
    maxNumberLostLayers = cms.uint32(2),
    minHitsToBypassChecks = cms.uint32(20),

    # Absolute cuts in case of no PV. If yes, please define also max_d0NoPV and max_z0NoPV 
    applyAbsCutsIfNoPV = cms.bool(False),
    max_d0NoPV = cms.double( 100.0 ),
    max_z0NoPV = cms.double( 100.0 ),

    # parameters for cutting on pterror/pt and number of valid hits
    max_relpterr = cms.double(9999.),
    min_nhits = cms.uint32(0),

    max_minMissHitOutOrIn = cms.int32(99),
    max_lostHitFraction = cms.double(1.0),

    # parameters for cutting on eta
    max_eta = cms.double(9999.),
    min_eta = cms.double(-9999.)

)


#process.elGsfTracksWithQuality.src = cms.InputTag("electronGsfTracks")

# sequence for re-running gsfTracking over RECO
process.myGsfReco = cms.Sequence(
    process.siPixelRecHits
    *process.siStripMatchedRecHits #make local hits
    *process.MeasurementTrackerEvent #ADD
    *process.siPixelClusterShapeCache # needed to add when moving from CMSSW_7_1_0_pre5 to pre7
    *process.iterTracking
    
    *process.electronSeedsSeq #ADD
    *process.electronSeeds    #produced merged collection of TkDriven and Ecaldriven seeds
    *process.electronCkfTrackCandidates
    *process.electronGsfTracks #run electron tracking
    *process.elGsfTracksWithQuality
    *process.printEventContent
)

outdir = "out_tests/"
#outfilename = outdir + "reGsfTracking_maxCand_" + str(maxCand) + "_MaxChi2_" + str(maxChi2) + "_nSigma_" + str(nSigma) + ".root"
outfilename = "trackValTree_reTrk.root"
print "Writing output to file: " + outfilename

process.TFileService = cms.Service("TFileService", # if save
                                   fileName = cms.string(outfilename)
                                   )



process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string("reco") #Quality_SimToReco = shared hits/#reco(or #sim)
process.TrackAssociatorByHits.Quality_SimToReco = cms.double(0.75)
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

process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.preValidation = cms.Sequence(
    process.simHitTPAssocProducer
    )

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

# paths
process.p = cms.Path(
    process.myGsfReco 
    *process.ValidationSelectors

    *process.preValidation
#    *process.printEventContent 
    *process.trackValTreeMaker
    )


process.schedule = cms.Schedule(
    process.p
    )
