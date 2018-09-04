#!/bin/csh -f
#
# Author: Greg Alushin galushin@rockefeller.edu
#
# This shell script performs post-processing on oversymmetrized microtubule reconstructions generated with 'HP' symmetry in frealign v9.11
#
# It regenerates microtubules from the correct protofilament of full and half-maps, then calculates FSC between half-maps
#
# Adapated from Liz Kellogg's mode0_script.com which interfaced with previous version of frealign (v9.09)
#
# Utilizes python scripts mtSeamSym.smart.py and mtFSC.py
#
# Assumes half-maps are in "scratch" directory, as is frealign default

set start 	= `grep start_process mparameters | awk '{print $2}'`
set data_input 	= `grep data_input mparameters | awk '{print $2}'`
set apix	= `grep pix_size mparameters | awk '{print $2}'`
set twist       = `grep ALPHA mparameters | awk '{print $2}'`
set rise        = `grep RISE mparameters | awk '{print $2}'`


mv ${data_input}_${start}_r1.mrc ${data_input}_${start}_r1_oversym.mrc

./mtSeamSym.smart.py -v ${data_input}_${start}_r1_oversym.mrc  --apix=${apix} --rise=${rise} --twist=${twist}
mv output.mrc ${data_input}_${start}_r1.mrc

mv scratch/${data_input}_${start}_r1_map1.mrc ${data_input}_${start}_r1_oversym_map1.mrc
./mtSeamSym.smart.py -v ${data_input}_${start}_r1_oversym_map1.mrc --apix=${apix} --rise=${rise} --twist=${twist}
mv output.mrc ${data_input}_${start}_r1_map1.mrc

mv scratch/${data_input}_${start}_r1_map2.mrc ${data_input}_${start}_r1_oversym_map2.mrc
./mtSeamSym.smart.py -v ${data_input}_${start}_r1_oversym_map2.mrc --apix=${apix} --rise=${rise} --twist=${twist}
mv output.mrc ${data_input}_${start}_r1_map2.mrc

./mtFSC.py -e ${data_input}_${start}_r1_map1.mrc -o ${data_input}_${start}_r1_map2.mrc --apix ${apix} --orad 145

echo "Post-processing of HP cycle ${start} finished.... "`date`
