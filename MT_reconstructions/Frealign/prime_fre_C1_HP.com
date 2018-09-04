#!/bin/bash

#this script will prime frealign by first calculating a C1 and an HP reconstruction from provided mparameters and stack
# It will then doing a single round of C1 refinement / C1 + HP reconstruction with parameter mask "0", priming scores without changing alignment parameters

#calculate initial reconstruction in C1
./editMparameters.py --msk 0 --num 1 --asym 0
frealign_calc_reconstructions 1
monitor_frealign.com
./postprocessing_script_C1.com

#calculate initial reconstruction in HP
./editMparameters.py --msk 0 --num 1 --asym HP
frealign_calc_reconstructions 1
monitor_frealign.com
./postprocessing_script_HP.com

#run one round of refinement with mask 0 in C1
./editMparameters.py --msk 0 --num 2 --asym C1
frealign_run_refine
monitor_frealign.com
./postprocessing_script_C1.com

#calculate HP reconstruction and correct seam from these alignment parameters
./editMparameters.py --msk 0 --num 2 --asym HP
frealign_calc_reconstructions 2
monitor_frealign.com
./postprocessing_script_HP.com
