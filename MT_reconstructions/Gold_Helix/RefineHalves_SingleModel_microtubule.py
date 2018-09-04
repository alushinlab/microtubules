#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will perform gold-standard FSC (i.e. independent half dataset) refinement 
# of a single model from a multimodel refinement
# Note that it is important to return to a highly-filtered initial model!!!!
#

#
#import modules, this could probably stand to be pruned
#

from __future__ import division

# import system / house keeping modules
import os
from optparse import OptionParser
import sys
import optparse
import glob
import subprocess
import linecache
import struct
import shutil
import time

# import eman / sparx specific modules for dealing with images
from EMAN2 import *
from sparx import *
from global_def import *

#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog --half_a <phase_flipped_half_stack_1> --half_b <phase_flipped_half_stack2> --init <initial_model> --apix <apix> --ang <angular steps> --dtrans <translational steps> --trans <translational range> --hpars <twist rise> --pfs <number of protofilaments> --nproc <num processors>")

	parser.add_option("--half_a",dest="stack_a",type="string",metavar="FILE",
		help="half stack 1 (output from HalfStacksfromMulti.py)")
	parser.add_option("--half_b",dest="stack_b",type="string",metavar="FILE",
		help="half stack 2 (output from HalfStacksfromMulti.py)")
	parser.add_option("--init",dest="init_mod",type="string",metavar="FILE",
		help="initial model (highly filtered)")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",default=2.18,
                help="angstroms per pixel")
	parser.add_option("--ang",dest="angular_steps",type="string",metavar="STRING",default='-1 -1 30 30 10 10',
		help="angular step series, separated by spaces and surrounded by single quotes")
	parser.add_option("--delta",dest="delta",type="string",metavar="STRING",default='4 4 3 3 2 2',
		help="angular step of reference projections, separated by spaces and surrounded by single quotes")
	parser.add_option("--dtrans",dest="trans_steps",type="string",metavar="STRING",default='4 4 2 2 1 1',
		help="translational step series, separated by spaces and surrounded by single quotes")
	parser.add_option("--trans",dest="trans_ranges",type="string",metavar="STRING",default='16 16 8 8 4 4',
		help="translational range series, separated by spaces and surrounded by single quotes")
	parser.add_option("--hpars",dest="helical_params",type="string",metavar="STRING",default='-27.66 9.88',
		help="twist and rise, separated spaces and surrounded by single quotes. Defaults for 13pf mts. Left handed twist is negative.")
	parser.add_option("--pfs",dest="pfs",type="int",metavar="int",default=13,
		help="Number of protofilaments.  Default 13")
	parser.add_option("--nprocs",dest="n_procs",type="int",metavar="INT",default=32,
		help="number of processors.  Should be an even integer multiple of 36 to run correctly on Rockefeller nodes.")
	parser.add_option("--crit_pass",dest="crit_pass",type="float",metavar="FLOAT",default=0.5,
		help="criterion for choosing butterworth lowpass filter pass frequency from FSC, default 0.5")
	parser.add_option("--crit_stop",dest="crit_stop",type="float",metavar="FLOAT",default=0.3,
		help="criterion for choosing butterworth lowpass filter stop frequency from FSC, default 0.3")
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params
#=========================
def checkFileExistence(params):
	if not os.path.isfile(params['stack_a']):
		print "\nError: stack file '%s' does not exist\n" % params['stack_a']
		sys.exit()
	if not os.path.isfile(params['stack_b']):
		print "\nError: stack file '%s' does not exist\n" % params['stack_b']
		sys.exit()
	if not os.path.isfile(params['init_mod']):
		print "\nError: 3D model '%s' does not exist\n" % params['init_mod']
		sys.exit()
#======================Microtubule-specific symmetry functions=======================
#===========================
def createWedgeMask(nx,rise,twist,apix,ovlp,wcmask=None):
	"""
	a soft wedge that follows helical symmetry
	"""
	import math
	img = EMData(nx,nx)
	img.to_zero()

	# find csym number from rotation
	csym=int(round(360.0/abs(twist)))

	#add ovlp degrees to overlap with the neighboring density
	overlap=ovlp*math.pi/180.0
	alpha = math.pi/2 - math.pi/csym
	for x,y in ((x,y) for x in range(0,int(nx)) for y in range(int(nx/2),int(nx))):
		dx = abs(x-nx/2)
		dy = abs(y-nx/2)
		# if above the line y = tan(alpha)*x
		inner = dx*math.tan(alpha)
		outer = dx*math.tan(alpha-overlap)
		if dy >= inner:
			img.set(x,y,1)
		elif dy >= outer:
			pos = (inner-dy)/(inner-outer)
			img.set(x,y,1-pos)

	img.process_inplace("mask.sharp",{"outer_radius":nx/2})

	wedge = EMData(nx,nx,nx)
	alpha = 360+(csym*twist)
	lrise = csym*rise
	rot = alpha/lrise*apix
	for z in range(nx):
		finalrot = ((z-nx/2)*rot)/3
		t = Transform()
		t.set_rotation({"type":"2d","alpha":-finalrot})
		newslice=img.process("xform",{"transform":t})
		wedge.insert_clip(newslice,(0,0,z))

	if wcmask is not None:
		# for additional masking of features inside wedge
		xmsk = int(wcmask[0]/apix)
		ymsk = int(wcmask[1]/apix)
		mskrad = int(wcmask[2]/apix)

		# see if mask is near the edge:
		edge=ymsk*math.atan(math.pi/csym)
		if (abs(xmsk)+mskrad)>=edge:
			# distance for corresponding positive mask
			edge = int(2*edge)
			xmsk2 = int(math.copysign(edge-abs(xmsk),xmsk)*-1)
			# take max of 1 mask
			avgr = Averagers.get("minmax",{"max":1})
			avgr.add_image_list([wedge,wedgeCylMask(nx,mskrad,xmsk2,ymsk,rot,pos=True)])
			wedge=avgr.finish()
		# multiply 0 mask
		wedge *= wedgeCylMask(nx,mskrad,xmsk,ymsk,rot)

	# odd-numbered protofilaments are off by 1/2 twist
	if csym%2==1:
		t = Transform({"type":"spider","psi":twist/2})
		wedge.process_inplace("xform",{"transform":t})

	#wedge.write_image('wedge_mask_p%d.mrc'%csym)

	return wedge
#===========================
def applySeamSym(vol,rise,rot,apix):
	"""
	apply seam symmetry based on results from Egelman search
	"""
	# find protofilament number from rotation
	sym=int(round(360.0/abs(rot)))

	rise/=apix
	# apply protofilament symmetry
	sumvol = vol.copy()
	pfoffset=int(sym/2)
	for pnum in range(-pfoffset,sym-pfoffset):
		if pnum==0: continue
		ang = rot*pnum
		trans = -(rise*pnum)
		t = Transform({"type":"spider","psi":ang})
		t.set_trans(0,0,trans)
		volcopy = vol.process("xform",{"transform":t})
		sumvol.add(volcopy)

	sumvol.process_inplace("normalize")
	return sumvol
#===========================
def regenerateFromPF(vol,wedgemask,rise,rot,apix):
	"""
	mask out one protofilament and regenerate the full microtubule
	"""
	#from reconstruction_rjh import smart_add
	# convert rise to pixels
	nx = vol.get_xsize()
	rise/=apix
	sym=int(round(360.0/abs(rot)))

	# apply protofilament symmetry
	sumvol = vol*wedgemask
	# save a copy of the single pf
	sumvol.write_image("pf.hdf",0)

	pfoffset=int(sym/2)
	for pnum in range(-pfoffset,sym-pfoffset):
		if pnum==0:
			continue
		ang = -(rot*pnum)
		trans = rise*pnum
		#print pnum, ang, trans
		t = Transform({"type":"spider","psi":ang})	
		t.set_trans(0,0,trans)
		volcopy = vol.process("xform",{"transform":t})
		seammaskcopy = wedgemask.process("xform",{"transform":t})
		sumvol = sumvol*(1-seammaskcopy)+volcopy*seammaskcopy

	sumvol.process_inplace("normalize")
	return sumvol
#================================================end microtubule-specific symmetry functions======================================


#=========================
def checkNprocSanity(params):
	num_procs=float(params['n_procs'])
	num_nodes = num_procs / 36
	integer_test = num_nodes.is_integer()

	if integer_test is False:
		print "nproc should be divisible by 36 to make best use of the cluster at Rockefeller.\n"
		sys.exit()

#=========================
def genClusterScript(runname,stackname,refname,nproc,outdir,xrad,t_s,delta,an,twist,rise,pfs,scriptout):
	#this module generates a run.sh style script which can be submitted to the Rockefeller cluster via SLURM
	shebang="""#!/bin/bash"""

	pbscommand1="""#SBATCH --job-name='"""
	pbscommand2="""'"""
	pbscommand3="""#SBATCH --export=ALL"""
	pbscommand4="""cd $SLURM_SUBMIT_DIR"""	
	pbscommand5="""PATH_TO_REFINE_SCRIPT='/ru-auth/local/home/sti/for_jove/scripts/EMAN2/EMAN2_actin2/refine.py'"""
	#pbscommand6="""source /usr/local/EMAN1-1.9/eman.bashrc"""
	#pbscommand7="""module load EMAN2/2.1"""

	executecommand1="""mpirun.eman2 -np """
	executecommand2=""" $PATH_TO_REFINE_SCRIPT """
	executecommand3=""" """
	executecommand4=""" """
	executecommand5=""" --ou=120 --rs=1 --xr='"""
	executecommand6="""' --ts='"""
	executecommand7="""' --delta='"""
	executecommand8="""' --an='"""
	executecommand9="""' --snr=0.08 --maxit=1 --ref_a=S --cutoff=10 --pix_cutoff=10000 --MPI --olmask=70 --ilmask=75 --lmask=280 --hpars='"""
	executecommand10=""" """
	executecommand11="""' --protos='"""
	executecommand12="""' --hsearch='90.0 180.0' --oplane=10 --recon_pad=2 --findseam --full_output > out.log"""
	
	
	f=open(scriptout, 'w')
	f.write("%s\n%s%s%s\n\n%s\n\n%s\n\n%s\n\n"%(shebang,pbscommand1,runname,pbscommand2,pbscommand3,pbscommand4,pbscommand5))
	f.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n"%(executecommand1,str(nproc),executecommand2,stackname,executecommand3,refname,executecommand4,outdir,executecommand5,str(xrad),executecommand6,str(t_s),executecommand7,str(delta),executecommand8,str(an),executecommand9,str(twist),executecommand10,str(rise),executecommand11,str(pfs),executecommand12))
	f.close()
#===========================
def SubmitandWait(script_tosubmit,nproc,delta_step,ref_dir):
	num_nodes = int(nproc)
	
	cmd = "sbatch --ntasks=%s --exclusive %s"%(str(num_nodes),script_tosubmit)
	print(cmd)
	subprocess.Popen(cmd,shell=True).wait()

	#get path to parameters output, which is last file written by refinement
	working_dir = os.getcwd()
	param_path="%s/%s/paramout_00%s_00"%(working_dir,ref_dir,str(delta_step))	
	
	#print param_path

	while not os.path.exists (param_path):
		time.sleep(10)
#============================
def GenHParsFile(twist,rise):
	hpars_file = "helical_params.spi"
	f = open(hpars_file,'w')
	f.write("; spi/spi\n")
	f.write("    1 2 %11.5f %11.5f\n"%(twist,rise))
	f.close()
#============================
def FindFSCcut(fsc_curve,cutoff):
	for i in range(len(fsc_curve[1])):
		if fsc_curve[1][i]<=cutoff:
			band=fsc_curve[0][i]
			break
	try:
		return (band)
	except NameError:
		print "Warning: FSC curve did not fall below %s cutoff."%(str(cutoff))
		return (False)
#============================
def PostRecon(vol_a,vol_b,mask,hpars,apix,recon_round,pass_crit,stop_crit):
	#sum vols for calculating helical parameters
	cmd_sum = "/programs/x86_64-linux/eman/1.9/bin/proc3d %s %s sum=tempsum.mrc"%(vol_a,vol_b)
	subprocess.Popen(cmd_sum,shell=True).wait()
	
	#
	# calculate new helical parameters
	#

	#rotate temp volume and convert to spider
	cmd_rot = "/programs/x86_64-linux/eman/1.9/bin/proc3d tempsum.mrc tempsumrot.spi rot=90,0,0 spidersingle"
	subprocess.Popen(cmd_rot,shell=True).wait()

	#run hsearch_lorentz
	cmd_hsearch = "hsearch_lorentz tempsumrot.spi %s %.4f 90.0 180.0 0.03 0.03"%(hpars,apix)
	subprocess.Popen(cmd_hsearch,shell=True).wait()
	
	#
	# impose helical parameters on each half volume
	#
	
	#get most recent rise and twist: code taken from Gabe Lander's helical module
	f = open(hpars)
	lines = f.readlines()
	f.close()
	last_pars=lines[-1].strip().split()
	
	#get rise and twist
	rise=float(last_pars[-1])
	twist=float(last_pars[-2])
	


	#open and helicize volumes
	
	a=EMData()
	a.read_image(vol_a)
	a=applySeamSym(a,rise,twist,apix)

	b=EMData()
	b.read_image(vol_b)
	b=applySeamSym(b,rise,twist,apix)

	a.write_image("%s_OverSym.hdf"%vol_a[0:-4])
	b.write_image("%s_OverSym.hdf"%vol_b[0:-4])	
	
	#make a wedge without internal features
	wcmask=None
	wedgemask=EMData()
	wedgemask=createWedgeMask(a.get_xsize(),rise,twist,apix,3,wcmask)

	#regenerate seamed microtubules
	a = regenerateFromPF(a,wedgemask,rise,twist,apix)
	a.write_image("%s_h.hdf"%vol_a[0:-4])
	b = regenerateFromPF(b,wedgemask,rise,twist,apix)
	b.write_image("%s_h.hdf"%vol_a[0:-4])

	#calculate FSC
	fsc_filename=("gold_fsc_round_%i"%(recon_round))
	
	#filter half volumes to use as reference for next round
	maskfile=EMData()
	maskfile.read_image(mask)
	w=1
	fsc_gold=fsc_mask(a,b,maskfile,w,fsc_filename)
	
	filt_pass=FindFSCcut(fsc_gold,pass_crit)
	filt_stop=FindFSCcut(fsc_gold,stop_crit)

	#also process sum volume for looking at

        fsc_point143=FindFSCcut(fsc_gold,0.143)
        filt143=float(apix/fsc_point143)

        print "FSC 0.143 between half volumes is %4.2f"%filt143
        mvcommand="mv tempsum.mrc sumvol_unfilt_round_%i.mrc"%recon_round
        subprocess.Popen(mvcommand,shell=True).wait()

	#go with butterworth filter
	try:	
		a_filtered=filt_btwl(a,filt_pass,filt_stop,True)
		b_filtered=filt_btwl(b,filt_pass,filt_stop,True)
		a_filtered.write_image("%s_h_filt.hdf"%vol_a[0:-4])
		b_filtered.write_image("%s_h_filt.hdf"%vol_b[0:-4])
		print "Reference half volumes filtered at pass band: %4.2f stop band: %4.2f"%((apix/filt_pass),(apix/filt_stop))
		return 
	except:
		#still write unfiltered references if finding cutoffs fails
		print "Warning!!!! no filter applied, continuing with misnamed unfiltered references"
		a_filtered=a
		b_filtered=b
		a_filtered.write_image("%s_h_filt.hdf"%vol_a[0:-4])
       		b_filtered.write_image("%s_h_filt.hdf"%vol_b[0:-4])
		return	
#
#=========================
#
if __name__ == "__main__":

	#initial setup and sanity check
	params=setupParserOptions()
	#loadEMANPaths()
	checkFileExistence(params)
	checkNprocSanity(params)

	#unpack some parameters to make things simpler
	init_mod=str(params['init_mod'])
	apix=float(params['apix'])
	nprocs=int(params['n_procs'])
	init_model=str(params['init_mod'])
	crit_pass=float(params['crit_pass'])
	crit_stop=float(params['crit_stop'])
	pfs=int(params['pfs'])

	angular_steps=str(params['angular_steps']).strip().split()
	trans_steps=str(params['trans_steps']).strip().split()
	trans_ranges=str(params['trans_ranges']).strip().split()
	deltas=str(params['delta']).strip().split()
	helical_params=str(params['helical_params']).strip().split()
	
	twist=float(helical_params[0])
	rise=float(helical_params[1])

	#generate helical parameters file
	GenHParsFile(twist,rise)
        hpars="helical_params.spi"
	
	#main loop
	
	for recon_round in range(len(angular_steps)):
		print "Beginning round %i"%recon_round
		#for simplicity's sake will currently run jobs in sequence rather than parallel

		halves=('a','b')
		
		#for the hell of it get current hpars to see if half vols give back the same helical parameters
		#once again use Gabe Lander's code
		f=open(hpars)
		lines=f.readlines()
		f.close()
		current_hpars=lines[-1].strip().split()
		
		for half in range(len(halves)):
			if recon_round == 0:
				current_ref = init_mod
			else:
				current_ref = "model_%s_%02d_h_filt.hdf"%(halves[half],(recon_round-1))	
			
			current_runname='run_%s_%02d'%(halves[half],recon_round)
			current_scriptname='run_%s_%02d.sh'%(halves[half],recon_round)
			current_outdir='out_%s_%02d'%(halves[half],recon_round)
			current_stack=str(params['stack_%s'%(halves[half])])
                        
			#TESTING123TESTING123
			#for now make fake directory and files
			#fake_dir_cmd = "mkdir %s"%(current_outdir)
			#print fake_dir_cmd
			#subprocess.Popen(fake_dir_cmd,shell=True).wait()
			#relevant_output = "volNoSym_%03d_00.hdf"%(int(deltas[round]))
			#fake_vol_cmd = "touch %s/%s"%(current_outdir,relevant_output)
			#subprocess.Popen(fake_vol_cmd,shell=True).wait()
			###############################################

			#hpars 2 should be twist, hpars 3 should be rise
			
			genClusterScript(current_runname,current_stack,current_ref,nprocs,current_outdir,trans_ranges[recon_round],trans_steps[recon_round],deltas[recon_round],angular_steps[recon_round],float(current_hpars[2]),float(current_hpars[3]),pfs,current_scriptname)
			print "Submitting refinement for half %s"%halves[half]

			SubmitandWait(current_scriptname,nprocs,deltas[recon_round],current_outdir)
			
			#link unsymmetrized output volume back to main directory for post-processing
			relevant_output = "volNoSym_%03d_00.hdf"%int(deltas[recon_round])
			link_cmd = "ln -s %s/%s model_%s_%02d.hdf"%(current_outdir,relevant_output,halves[half],recon_round)
			subprocess.Popen(link_cmd,shell=True).wait()
			
			#link mask back to main directory, need not do every time but won't hurt
			link_mask_cmd = "ln -s %s/mask3D_cyl.mrc ."%(current_outdir)
			subprocess.Popen(link_mask_cmd,shell=True).wait()
	
		#post processing
		volume_a="model_a_%02d.hdf"%recon_round
		volume_b="model_b_%02d.hdf"%recon_round
		mask="mask3D_cyl.mrc"

		PostRecon(volume_a,volume_b,mask,hpars,apix,recon_round,crit_pass,crit_stop)
