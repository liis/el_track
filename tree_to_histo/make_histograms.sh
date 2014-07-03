#!/bin/bash

INFILE=$1 #file of the input tree to process

echo $INFILE
python make_eff_fake_hists.py --infile=$INFILE


