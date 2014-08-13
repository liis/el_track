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
