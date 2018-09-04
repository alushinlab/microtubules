#!/bin/bash

# This script will run an interative frealign C1/HP refinement between the specified rounds
# 
# first command-line argument is first round
# second command-line argument is last round
#
# Note that both a reconstruction and parameter file from the previous round MUST EXIST
# These should be generated with prime_fre_C1_HP.com if starting a new Frealign refinement from EMAN2

start=$1
end=$2
end=$((end+1))

while [ $start -lt $end ]
do
	./editMparameters.py --num $start --asym 0
	frealign_run_refine
	monitor_frealign.com
	./postprocessing_script_C1.com

	./editMparameters.py --num $start --asym HP
	frealign_calc_reconstructions $start
	monitor_frealign.com
	./postprocessing_script_HP.com
	
	start=$((start+1))
done
