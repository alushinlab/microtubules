#!/bin/csh
#unlimit
limit coredumpsize 0
set bindir = /home/galushin/frealign_v9.07/bin
set working_directory = `pwd`
set SCRATCH = ../scratch
setenv NCPUS 12

#
set start = $3
@ prev = $start - 1
set data_input	= $4
set raw_images	= $5
set thresh	= `grep thresh_reconst mparameters_run | awk '{print $2}'`
set pbc		= `grep PBC mparameters_run | awk '{print $2}'`
set boff	= `grep BOFF mparameters_run | awk '{print $2}'`
set dang	= `grep DANG mparameters_run | awk '{print $2}'`
set itmax	= `grep ITMAX mparameters_run | awk '{print $2}'`
set mode	= `grep MODE mparameters_run | awk '{print $2}'`
set FFILT	= `grep FFILT mparameters_run | awk '{print $2}'`
set FBEAUT	= `grep FBEAUT mparameters_run | awk '{print $2}'`
set rrec	= `grep res_reconstruction mparameters_run | awk '{print $2}'`
set rref	= `grep res_high_refinement mparameters_run | awk '{print $2}'`
set rbf		= `grep RBfactor mparameters_run | awk '{print $2}'`
set sym		= `grep Sym mparameters_run | awk '{print $2}'`
set alpha	= `grep ALPHA mparameters_run | awk '{print $2}'`
set rise	= `grep RISE mparameters_run | awk '{print $2}'`
set nsub	= `grep NSUBUNITS mparameters_run | awk '{print $2}'`
set nstart	= `grep NSTARTS mparameters_run | awk '{print $2}'`
set stiff	= `grep STIFFNESS mparameters_run | awk '{print $2}'`
set pix		= `grep pix_size mparameters_run | awk '{print $2}'`
set kV		= `grep kV mparameters_run | awk '{print $2}'`
set AmpC	= `grep Amp_contrast mparameters_run | awk '{print $2}'`
set dstep	= `grep dstep mparameters_run | awk '{print $2}'`
set ro		= `grep outer_radius mparameters_run | awk '{print $2}'`
set ri		= `grep inner_radius mparameters_run | awk '{print $2}'`
set MW		= `grep mol_mass mparameters_run | awk '{print $2}'`
set cs		= `grep Cs mparameters_run | awk '{print $2}'`

if ($mode == 0) then
  cp ${data_input}_${prev}.mrc ${data_input}_${start}.mrc
  set parm = ${data_input}_${prev}.par

else

  set parm = ${data_input}_${start}.par

endif
#
cd $SCRATCH
#
\rm ${data_input}_${start}.res
#
time ${bindir}/frealign_v9_mp.exe << eot >& ${data_input}_mreconstruct.log
M,0,F,F,F,F,0,${FBEAUT},${FFILT},F,F,0,F,2,1			!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FSTAT,IMEM,INTERP
${ro},${ri},${pix},${MW},${AmpC},0.0,${pbc},0.0,10.,1,10	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
1 1 1 1 1							!MASK
${1},${2}							!IFIRST,ILAST 
${sym}								!ASYM symmetry card (I=icosahedral)
${alpha},${rise},${nsub},${nstart},${stiff}			!ALPHA,RISE,NSUBUNITS,NSTARTS,STIFFNESS
1.,${dstep},60.0,${thresh},${cs},${kV},0.0,0.0			!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
${rrec}, 200.0, ${rref}, 100.0, ${rbf}				!RREC,RMIN,RMAX,DFSTD,RBFACT
${working_directory}/${raw_images}.mrc
/dev/null
${working_directory}/${parm}
${data_input}_${start}.res
${data_input}_${start}_dummy.shft
0., 0., 0., 0., 0., 0., 0., 0.					! terminator with RELMAG=0.0
${data_input}_${start}.mrc
${data_input}_${start}_weights
${data_input}_${start}_map1.mrc
${data_input}_${start}_map2.mrc
${data_input}_${start}_phasediffs
${data_input}_${start}_pointspread
eot
#
mv ${data_input}_${start}.mrc ${working_directory}/.
\rm ${data_input}_${start}_weights
\rm ${data_input}_${start}_dummy.shft
#\rm ${data_input}_${start}_map1.mrc
#\rm ${data_input}_${start}_map2.mrc
\rm ${data_input}_${start}_phasediffs
# \rm ${data_input}_${start}_pointspread
#
echo 'mreconstruct.com finished' >> ${data_input}_mreconstruct.log
date
#
