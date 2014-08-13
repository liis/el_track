ALL_SUBSYSTEMS+=$(patsubst src/%,%,src/MakeTree)
subdirs_src_MakeTree = src_MakeTree_MakeTrackValTree
ALL_PACKAGES += $(patsubst src/%,%,src/MakeTree/MakeTrackValTree)
subdirs_src_MakeTree_MakeTrackValTree := src_MakeTree_MakeTrackValTree_python src_MakeTree_MakeTrackValTree_src
ifeq ($(strip $(PyMakeTreeMakeTrackValTree)),)
PyMakeTreeMakeTrackValTree := self/src/MakeTree/MakeTrackValTree/python
src_MakeTree_MakeTrackValTree_python_parent := 
ALL_PYTHON_DIRS += $(patsubst src/%,%,src/MakeTree/MakeTrackValTree/python)
PyMakeTreeMakeTrackValTree_files := $(patsubst src/MakeTree/MakeTrackValTree/python/%,%,$(wildcard $(foreach dir,src/MakeTree/MakeTrackValTree/python ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
PyMakeTreeMakeTrackValTree_LOC_USE := self  
PyMakeTreeMakeTrackValTree_PACKAGE := self/src/MakeTree/MakeTrackValTree/python
ALL_PRODS += PyMakeTreeMakeTrackValTree
PyMakeTreeMakeTrackValTree_INIT_FUNC        += $$(eval $$(call PythonProduct,PyMakeTreeMakeTrackValTree,src/MakeTree/MakeTrackValTree/python,src_MakeTree_MakeTrackValTree_python,1,1,$(SCRAMSTORENAME_PYTHON),$(SCRAMSTORENAME_LIB),,))
else
$(eval $(call MultipleWarningMsg,PyMakeTreeMakeTrackValTree,src/MakeTree/MakeTrackValTree/python))
endif
ALL_COMMONRULES += src_MakeTree_MakeTrackValTree_python
src_MakeTree_MakeTrackValTree_python_INIT_FUNC += $$(eval $$(call CommonProductRules,src_MakeTree_MakeTrackValTree_python,src/MakeTree/MakeTrackValTree/python,PYTHON))
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
ALL_SUBSYSTEMS+=$(patsubst src/%,%,src/TestReReco)
subdirs_src_TestReReco = src_TestReReco_submission_scripts
ALL_PACKAGES += $(patsubst src/%,%,src/TestReReco/submission_scripts)
subdirs_src_TestReReco_submission_scripts := src_TestReReco_submission_scripts_templates
