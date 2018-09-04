#!/usr/bin/env python

import os,sys
import subprocess
import optparse
import math
from itertools import product

#==========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -s <stack>")
	parser.add_option("-e",dest="even",type="string",metavar="FILE",
		help="even volume")
	parser.add_option("-o",dest="odd",type="string",metavar="FILE",
		help="odd volume")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
		help="pixel size")
	parser.add_option("--orad", dest="orad", type="int", metavar="#", default=145,
		help="outer radius (in Angstroms, default=145)")
	parser.add_option("--irad", dest="irad", type="int", metavar="#", default=75,
		help="inner radius (in Angstroms, default=75)")
	parser.add_option("--lmask", dest="lmask", type="int", metavar="#", default=240,
		help = "length of microtubule to keep (in Angstroms, default=240)")

	options,args = parser.parse_args()

	if len(args) > 1:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()
	params={}

	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params

#==========================
def checkConflicts(params):
	if not params['even'] or not params['odd']:
		print "Specify an even and odd volume"
		sys.exit()
	if not os.path.isfile(params['even']):
		print "the specified volume '%s' does not exist"%params['even']
		sys.exit()
	if not os.path.isfile(params['odd']):
		print "the specified volume '%s' does not exist"%params['odd']
		sys.exit()
	if not params['apix']:
		print "Specify a pixel size"
		sys.exit()
	# read volume and get size
	print "reading file: %s"%params['even']
	params['evol'] = EMData()
	params['evol'].read_image(params['even'])
	params['nx'] = params['evol'].get_xsize()
	print "reading file: %s"%params['odd']
	params['ovol'] = EMData()
	params['ovol'].read_image(params['odd'])

# Generate 2D slices to be inserted into mask3D volume
def createMask2D(params):
	from itertools import product
	apix = params['apix']
	orad = float(params['orad'])/apix
	irad = float(params['irad'])/apix
	nx = params['nx']
	falloff_r = 20			# use steeper falloff
	mask2D = EMData(nx,nx)
	mask2D.to_one()
	for x,y in product(range(nx),range(nx)):
		dx = abs(x-nx/2)
		dy = abs(y-nx/2)
		r2 = dx**2+dy**2
		if r2 > orad*orad:
			wt1 = 0.5*(1 + math.cos(math.pi*min(1,(math.sqrt(r2)-orad)/falloff_r)))
			mask2D.set(x,y,wt1)
		elif r2 < irad*irad:
			wt2 = 0.5*(1 + math.cos(math.pi*min(1,(irad-math.sqrt(r2))/falloff_r)))
			mask2D.set(x,y,wt2)
	#mask2D.write_image('mask.mrc')
	return mask2D

def createMask3D(params,mask2D):
	apix = params['apix']
	orad = float(params['orad'])/apix
	irad = float(params['irad'])/apix
	nx = mask2D.get_xsize()
	
	zrad = float(params['lmask'])/apix/2
	
	mask3D = EMData(nx,nx,nx)
	falloff_z = 30.0
	# now apply soft mask
	for z in range(nx):
		img = EMData(nx,nx)
		img = mask2D.copy()
		# here "img = mask2D" won't work !!
		dz = abs(z-nx/2)
		if dz > zrad:
			wt3 = 0.5*(1 + math.cos(math.pi*min(1,(dz-zrad)/falloff_z)))
			img.mult(wt3)
			#img.write_image("test_%d.mrc"%z)
		mask3D.insert_clip(img,(0,0,z))
	
	#mask3D.write_image('mask3D_%dx%dx%d.mrc'%(orad,irad,zrad))
	return mask3D

#===========================
def createCylMask(params):
	"""
	create a cylindrical mask with gaussian edges
	"""
	print "creating mask"
	#from itertools import product
	import math

	apix = params['apix']
	nx = params['nx']

	## convert mask values to pixels
	rmax = int(params['orad']/apix)
	rmin = int(params['irad']/apix)
	lmask = int((params['lmask']/apix)/2)
	falloff_outer = lmask*0.4
	falloff_inner = rmin*0.4

	## first create cylinder with inner & outer mask
	cyl = EMData(nx,nx,nx)
	for i in range(nx):
		mask=EMData(nx,nx)
		mask.to_one()
		## mask the inner & outer radii
		#for x,y in ((x,y) for x in range(nx) for y in range(nx)):
		for x,y in product(range(nx),range(nx)):
			dx = abs(x-nx/2)
			dy = abs(y-nx/2)
			r2 = dx**2+dy**2
			if r2 > rmax*rmax:
				wt1 = 0.5*(1 + math.cos(math.pi*min(1,(math.sqrt(r2)-rmax)/falloff_outer)))
				mask.set(x,y,wt1)
			elif r2 < rmin*rmin:
				wt2 = 0.5*(1 + math.cos(math.pi*min(1,(rmin-math.sqrt(r2))/falloff_inner)))
				mask.set(x,y,wt2)
		## mask along length
		dz = abs(i-nx/2)
		if dz > lmask:
			wt3 = 0.5*(1+math.cos(math.pi*min(1,(dz-lmask)/falloff_outer)))
			mask.mult(wt3)
		cyl.insert_clip(mask,(0,0,i))
	
	cyl.write_image('cyl_fsc.mrc')
	return cyl

#==========================
def getEMANPath():
	### get the eman2 directory	
	emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
	if emanpath:
		emanpath = emanpath.replace("EMAN2DIR=","")
	if os.path.exists(emanpath):
		return emanpath
	print "EMAN2 was not found, make sure it is in your path"
	sys.exit()

#==========================
def calcFSC(params,cyl):
	print "calculating FSC"
	fscc = fsc_mask(params['evol'],params['ovol'],cyl,1,"%s.fsc_orad%d"%(params['odd'][:-9],params['orad']))
	getResFromFSCFile("%s.fsc_orad%d"%(params['odd'][:-9],params['orad']),params['apix'])
	getResFromFSCFile("%s.fsc_orad%d"%(params['odd'][:-9],params['orad']),params['apix'],cutoff=0.143)

#==========================
def getResFromFSCFile(fscfile,apix,cutoff=0.5):
	f = open(fscfile, 'r')
	f2 = open('%s.fscNumber_orad%d'%(params['odd'][:-9],params['orad']),"a")
	lastx=0
	lasty=0
	for line in f:
		xy = line.strip().split()
		x = float(xy[0])
		y = float(xy[1])
		if x < 0.0 or x > 0.5:
			print "FSC is wrong data format"
		if y > cutoff:
			print "y not greater than cutoff %.3f %.2f"%(y,cutoff)
			#store values for later
			lastx = x
			lasty = y
		else:
			# get difference of fsc
			diffy = lasty-y
			# get distance from 0.5
			distfsc = (cutoff-y) / diffy
			# get interpolated spatial freq
			intfsc = x - distfsc * (x-lastx)
			print "interpolated spatial frequency is: " + str(intfsc) + "\n"
			# convert to Angstroms
			res = apix / intfsc
			print "Resolution at FSC=%.3f: %.2f"%(cutoff,res)
			f2.write("Resolution at FSC=%.3f: %.2f\n"%(cutoff,res))	
			break
	f.close()
	f2.close()

#==========================
if __name__ == "__main__":
	#getEMANPath()
	from EMAN2 import *
	from sparx import *
	params=setupParserOptions()
	checkConflicts(params)
	#try:
	#	cyl = EMData('cyl_fsc.mrc')
	#except:
	#       print "cannot find existing cyl_fsc.mrc, will create a new one"
	#	cyl = createCylMask(params)
	#cyl = createCylMask(params)
	mask2D = createMask2D(params)
	maskfsc = createMask3D(params,mask2D)
	maskfsc.write_image('maskFSC.mrc')
	calcFSC(params,maskfsc)
