ifeq ($(strip $(MakeTree/MakeTrackValTree)),)
ALL_COMMONRULES += src_MakeTree_MakeTrackValTree_src
src_MakeTree_MakeTrackValTree_src_parent := MakeTree/MakeTrackValTree
src_MakeTree_MakeTrackValTree_src_INIT_FUNC := $$(eval $$(call CommonProductRules,src_MakeTree_MakeTrackValTree_src,src/MakeTree/MakeTrackValTree/src,LIBRARY))
MakeTreeMakeTrackValTree := self/MakeTree/MakeTrackValTree
MakeTree/MakeTrackValTree := MakeTreeMakeTrackValTree
MakeTreeMakeTrackValTree_files := $(patsubst src/MakeTree/MakeTrackValTree/src/%,%,$(wildcard $(foreach dir,src/MakeTree/MakeTrackValTree/src ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
MakeTreeMakeTrackValTree_BuildFile    := $(WORKINGDIR)/cache/bf/src/MakeTree/MakeTrackValTree/BuildFile
MakeTreeMakeTrackValTree_LOC_USE := self  FWCore/Framework FWCore/PluginManager FWCore/ParameterSet DataFormats/TrackReco DataFormats/SiPixelDetId CommonTools/UtilAlgos DQMServices/Core SimDataFormats/TrackerDigiSimLink DataFormats/SiStripDetId Geometry/TrackerGeometryBuilder Geometry/Records MagneticField/Records MagneticField/Engine SimDataFormats/Vertex SimDataFormats/TrackingAnalysis DataFormats/EgammaReco SimTracker/Records RecoLocalTracker/ClusterParameterEstimator rootcintex root SimTracker/TrackAssociation SimTracker/TrackerHitAssociation DQMServices/ClientConfig Geometry/CommonDetUnit
MakeTreeMakeTrackValTree_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,MakeTreeMakeTrackValTree,MakeTreeMakeTrackValTree,$(SCRAMSTORENAME_LIB),src/MakeTree/MakeTrackValTree/src))
MakeTreeMakeTrackValTree_PACKAGE := self/src/MakeTree/MakeTrackValTree/src
ALL_PRODS += MakeTreeMakeTrackValTree
MakeTreeMakeTrackValTree_INIT_FUNC        += $$(eval $$(call Library,MakeTreeMakeTrackValTree,src/MakeTree/MakeTrackValTree/src,src_MakeTree_MakeTrackValTree_src,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
endif
