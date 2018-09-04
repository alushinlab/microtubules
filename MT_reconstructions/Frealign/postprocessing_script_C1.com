#!/bin/csh -f
#
# Author: Greg Alushin galushin@rockefeller.edu
#
# This shell script moves around the results of C1 reconstructions from frealign v9.11
#
# Preventing them from being overwritten by HP reconstructions from the same round

set start 	= `grep start_process mparameters | awk '{print $2}'`
set data_input 	= `grep data_input mparameters | awk '{print $2}'`

mv ${data_input}_${start}_r1.mrc ${data_input}_${start}_r1_C1.mrc

grep C ${data_input}_${start}_r1.par > ${data_input}_${start}_r1_C1.fsc

echo "C1 reconstruction and FSC of cycle ${start} retained.... "`date`
