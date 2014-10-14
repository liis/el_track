#!/bin/bash                                                                                                                                                                              

#INFILE=$1 #file of the input tree to process                                                                                                                                             
#OUTDIR=$2

#echo $INFILE
python resolution_plots_paramScans.py --indir=../tree_to_histo/histograms/13TeV_011014/
#python make_eff_fake_hists.py --infile=$INFILE --outdir=$OUTDIR