#!/bin/csh 
#unlimit
limit coredumpsize 0
set working_directory = `pwd`
set SCRATCH = ../scratch

set bindir = /home/galushin/frealign_v9.07/bin

set start = $3
set data_input	= $4
set raw_images	= $5
set target	= `grep thresh_refine mparameters_run | awk '{print $2}'`
set pbc		= `grep PBC mparameters_run | awk '{print $2}'`
set boff	= `grep BOFF mparameters_run | awk '{print $2}'`
set dang	= `grep DANG mparameters_run | awk '{print $2}'`
set itmax	= `grep ITMAX mparameters_run | awk '{print $2}'`
set mode	= `grep MODE mparameters_run | awk '{print $2}'`
set FPART	= `grep FPART mparameters_run | awk '{print $2}'`
set FMATCH	= `grep FMATCH mparameters_run | awk '{print $2}'`
set rrec	= `grep res_reconstruction mparameters_run | awk '{print $2}'`
set rref	= `grep res_high_refinement mparameters_run | awk '{print $2}'`
set rlowref	= `grep res_low_refinement mparameters_run | awk '{print $2}'`
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
set XSTD	= `grep XSTD mparameters_run | awk '{print $2}'`
set dstep	= `grep dstep mparameters_run | awk '{print $2}'`
set ro		= `grep outer_radius mparameters_run | awk '{print $2}'`
set ri		= `grep inner_radius mparameters_run | awk '{print $2}'`
set MW		= `grep mol_mass mparameters_run | awk '{print $2}'`
set cs		= `grep Cs mparameters_run | awk '{print $2}'`
set mode	= `grep MODE mparameters_run | awk '{print $2}'`
set mask	= `grep MASK mparameters_run | awk '{print $2}'`

@ prev = $start - 1
cd $SCRATCH

cp ${working_directory}/${data_input}_${prev}.par ${data_input}_${prev}.par_$1_$2
cp ${working_directory}/${data_input}_${prev}.mrc ${data_input}_${start}.mrc_$1_$2
\rm ${data_input}.par_${1}_${2} >& /dev/null

#

time ${bindir}/frealign_v9.exe << eot >& ${data_input}_mrefine_n.log_${1}_${2}
M,${mode},F,F,F,${FPART},0,F,F,F,${FMATCH},0,F,0,0			!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FSTAT,IMEM,INTERP
${ro},${ri},${pix},${MW},${AmpC},${XSTD},${pbc},${boff},${dang},${itmax},20	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
${mask},${mask},${mask},${mask},${mask}					!MASK
${1},${2}								!IFIRST,ILAST
${sym}									!ASYM symmetry card (I=icosahedral)
${alpha},${rise},${nsub},${nstart},${stiff}				!ALPHA,RISE,NSUBUNITS,NSTARTS,STIFFNESS
1.0,${dstep},${target},0.0,${cs},${kV},0.0,0.0				!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
${rrec},   ${rlowref},   ${rref}, 100.0, ${rbf}				!RREC,RMIN,RMAX,DFSTD,RBFACT
${working_directory}/${raw_images}.mrc
${data_input}_match.mrc_${1}_${2}
${data_input}_${prev}.par_$1_$2
${data_input}_${start}.par_${1}_${2}
${data_input}_${start}.shft_${1}_${2}
-100., 0., 0., 0., 0., 0., 0., 0.					! terminator with RELMAG=0.0
${data_input}_${start}.mrc_${1}_${2}
${data_input}_weights_${start}_${1}_${2}
${data_input}_map1_${1}_${2}
${data_input}_map2_${1}_${2}
${data_input}_phasediffs_${1}_${2}
${data_input}_pointspread_${1}_${2}
eot
#
