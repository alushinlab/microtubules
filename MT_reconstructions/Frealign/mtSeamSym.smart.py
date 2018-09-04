#!/usr/bin/env python

import os,sys
import subprocess
import optparse
import math
import time
from itertools import product

#==========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -s <stack>")
	parser.add_option("-v",dest="volume",type="string",metavar="FILE",
		help="input volume for masking")
	parser.add_option("--applysym",dest="applysym",action="store_true",default=False,
		help="apply real-space helical symmetry")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
		help="pixel size")
	parser.add_option("--rise",dest="rise", type="float",metavar="FLOAT",
		help="hpar rise")
	parser.add_option("--twist",dest="twist", type="float",metavar="FLOAT",
		help="hpar twist")
	parser.add_option("--vrise",dest="vrise", type="float",metavar="FLOAT",
		help="vertical rise")
	parser.add_option("--vtwist",dest="vtwist", type="float",metavar="FLOAT",
		help="vertical twist")
	parser.add_option("--orad", dest="orad", type="int", metavar="INT",default=215,
		help="maximum radius in x-y plane, in angstrom")
	parser.add_option("--irad", dest="irad", type="int", metavar="INT", default=50,
		help="minimum radius in x-y plane, in angstrom")
	parser.add_option("--zrad", dest="zrad", type="int", metavar="INT",
		help="radius in z direction, in angstrom")
	parser.add_option("--fast",dest="fast",action="store_true",default=False,
		help="use fast masking for testing")
	parser.add_option("--opposite",dest="opposite",action="store_true",default=False,
		help="seam is opposite of frealign")
	parser.add_option("--decor",dest="decor",type="choice", metavar="['kinesin','none','EB']",
		choices=['kinesin','none','EB'], default = 'kinesin',
		help="bound with kinesin or EB or nothing")

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
	if not params['volume']:
		print "Specify a volume"
		sys.exit()
	if not os.path.isfile(params['volume']):
		print "the specified volume '%s' does not exist"%params['volume']
		sys.exit()
	# read volume and get size
	print "reading file: %s"%params['volume']
	params['vol'] = EMData()
	params['vol'].read_image(params['volume'])
	params['nx'] = params['vol'].get_xsize()

	# get helical parameters
	params['pf']=int(round(360.0/abs(params['twist'])))
        params['finalrot'] = 1000

def cart2pol(x, y):
        import numpy as np
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return(rho, phi)

def pol2cart(th,r):
        import numpy as np
        x = r*np.cos(th*np.pi/180)
        y = r*np.sin(th*np.pi/180)
        return(x,y)

def sumdiff(in1,in2):
        sdiff = 0
        for i in range(len(in1)):
                sdiff += (in1[i]-in2[i])
        return sdiff


def initialize_assignments(n_elem,ncluster,mat):
        import numpy as numpy
        #assn = numpy.floor(numpy.random.rand(n_elem)*int(ncluster))
        assn = numpy.zeros(n_elem)
        theta = float(360)/float(ncluster)
        th = 0
        rad = 145
        pseudo_centers = numpy.zeros(shape=(ncluster,2))
        ii = 0
        c = float(params['nx'])/float(2)
        while ii < ncluster:
                pseudo_centers[ii,:] = pol2cart(th,rad)
                pseudo_centers[ii,:] += (c,c)
                ii+=1
                th += theta
        assn = assign_to_cluster(mat,pseudo_centers)
        return assn

def compute_centers(mat,assn,nexpect):
        import numpy as numpy
        n_elem = len(mat[:,0])
        n_cluster = len(numpy.unique(assn))
        cluster_centers = numpy.zeros(shape=(nexpect,2))
        cluster_sizes = numpy.zeros(shape=(nexpect,1))
        if(len(mat[:,0]) != len(assn)):
                print "ERROR! matrix should be the same size as assignment array, but they are not! " + str(len(mat[:,0])) + " " + str(len(assn)) + "\n"
                
        for j in range(nexpect):
                j_cluster = assn == j #find all indices which belong to this cluster
                cluster_sz = 0
                for jj in range(len(j_cluster)):
                        if j_cluster[jj]:
                                cluster_centers[j,0] += mat[jj,0]
                                cluster_centers[j,1] += mat[jj,1]
                                cluster_sz += 1
                cluster_sizes[j] = cluster_sz
                if cluster_sz >= 1:
                        cluster_centers[j,0] = cluster_centers[j,0] / cluster_sz
                        cluster_centers[j,1] = cluster_centers[j,1] / cluster_sz
                else:
                        cluster_centers[j,0] = 1000
                        cluster_centers[j,1] = 1000


        #this function doesn't fucking work!
        #sorted_size_ndx = cluster_sizes.argsort()
        max_sz = 0
        max_ndx = 0
        for ii in range(len(cluster_sizes)):
                if cluster_sizes[ii] > max_sz:
                        max_sz = cluster_sizes[ii]
                        max_ndx = ii
        
        for i in range(n_cluster):
                if cluster_sizes[i] == 0:
                        ndx_x = find_all_ndc(assn,max_ndx)
                        t = int(floor(numpy.random.rand(1)*len(ndx_x)))
                        #print "random index: " + str(t) + " len " + str(len(ndx_x)) + "\n"
                        randndx = ndx_x[t]
                        assn[randndx] = i
                        cluster_centers[i,0] = mat[randndx,0]
                        cluster_centers[i,1] = mat[randndx,1]
                else:
                        break
                        
        return cluster_centers
                
        
def assign_to_cluster(mat,centr):
        import numpy as numpy
        n_clusters = len(centr)
        n_elem = int(len(mat[:,0]))
        assn = numpy.zeros(shape=(n_elem,1))
        for ii in range(n_elem):
                best_cluster = 0
                best_dist = 10000
                for jj in range(n_clusters):
                        dd = numpy.linalg.norm( centr[jj,:] - mat[ii,:] )
                        if dd < best_dist:
                                best_dist = dd
                                best_cluster = jj
                assn[ii] = best_cluster
        return assn

def mean_distance(mat,assn,cntr):
        import numpy as numpy
        meand = numpy.zeros(shape=(len(cntr[:,0]),1))
        sz = numpy.zeros(shape=(len(cntr[:,0]),1))
        #print "number of unique assignments " + len(numpy.unique(assn)) + "\n"
        for ii in range(len(mat[:,0])):
                #print "ndx " + str(ii) + " has assn " + str(int(assn[ii])) +  "\n"
                meand[int(assn[ii])] += numpy.linalg.norm(mat[ii,:]-cntr[int(assn[ii]),:])
                sz[int(assn[ii])] += 1
        for jj in range(len(cntr[:,])):
                if sz[jj] > 0:
                        meand[jj] = meand[jj]/sz[jj]
                else:
                        print "how is this possible.. for " + str(jj) + " size is 0?? " + str(sz[jj]) + "\n"
        return meand
                

def kmeans_cluster(mat,nclusters):
        import numpy as numpy
        #randomly assign indices
        n_elem = int(len(mat[:,0]))
        curr_assn = initialize_assignments(n_elem,nclusters,mat)
        prev_assn = numpy.zeros(n_elem)
        #while cluster assignments still change
        itr = 0
        while(abs(sumdiff(curr_assn,prev_assn)) != 0 ):
                itr+=1
                print "on iteration " + str(itr) + "\n"
                #determine closest center
                prev_assn = curr_assn
                cluster_cntrs = compute_centers(mat,curr_assn,nclusters)
                #update assignments
                curr_assn = assign_to_cluster(mat,cluster_cntrs)                
        return curr_assn

#find the indices ii in assn
def find_all_ndc(assn,ii):
        import numpy as numpy
        ndx_i = []
        for indx in range(len(assn)):
                if assn[indx] == ii:
                        ndx_i = numpy.append(ndx_i,indx)
        return ndx_i

def locateGoodPF(params):
        import numpy as numpy

        if params['finalrot'] != 1000 or params['finalrot'] == -1:
                return params['finalrot']

        
        #first, copy the map and low-pass filter
        lpvol = params['vol'].process("filter.lowpass.gauss",{"cutoff_freq":0.1})
        #compute projection before determining threshold
        zproj_vol = lpvol.project("standard",Transform())
	zproj_vol.write_image("proj_vol.mrc")
        stdev = zproj_vol.get_attr("sigma")
        #print "std-dev is " + str(stdev) + "\n"
        px = zproj_vol.calc_highest_locations(3*stdev)
	if  len(px) == 0:
		print "[ERROR] there are no data-points in your volume! Are you sure this is a good reconstruction??\n"
		exit(1)
        xycoords = numpy.zeros(shape=(len(px),2))
        for pp in range(len(px)):
                pp_pt = px[pp].get_point()
                xycoords[pp,0] = pp_pt[0]
                xycoords[pp,1] = pp_pt[1]
        npf = params['pf']
        cluster_assn = kmeans_cluster(xycoords,npf)
        c = float(params['nx'])/float(2)
        cent = [c, c]

        #for ii in range(len(px)):
        #            print "assn: " + str(cluster_assn[ii]) + " " + str(ii) + " " + str(xycoords[ii,0]) + " " + str(xycoords[ii,1]) + "\n"

        
        pf_theta = numpy.zeros(npf)
        for  ii in range(npf):
		#print "clustering " + str(len(cluster_assn)) + " elements\n"
                ndx = find_all_ndc(cluster_assn,ii)
		#print "looking for index "+ str(ii) + " n-found: " + str(len(ndx)) + "\n"
                sum_theta = 0
                for jj in range(len(ndx)):
                    [theta,rad] = cart2pol(xycoords[jj,0]-cent[0],xycoords[jj,1]-cent[1])
                    sum_theta =+ (theta*180/numpy.pi)-90
                avg_theta = sum_theta / len(ndx)
                pf_theta[ii] = avg_theta
        
        #find theta closest to 0.
        best_theta = 180
        best_ndx = 0
        for ii in range(len(pf_theta)):
                #find closest positive angle 
                if pf_theta[ii] > 0 and pf_theta[ii] < best_theta :
                    best_theta = pf_theta[ii]
                    best_ndx = ii
        print "best-angle detected: " + str(best_theta) + "\n"
        params['finalrot'] = best_theta
        return best_theta
        
def circularMask2D(nx):
	falloff = 30.0					# cosine falloff on the edge
	smask2D = EMData(nx,nx)
	smask2D.to_one()
	rad = nx/2-falloff

	for i,j in product(range(nx),range(nx)):
			dx = abs(i-nx/2)
			dy = abs(j-nx/2)

			r2 = dx**2 + dy**2
			if r2 > rad**2:
				wt = 0.5*(1 + math.cos(math.pi*min(1,(math.sqrt(r2)-rad)/falloff)))
			else:
				wt = 1
			smask2D.set(i,j,wt)
	#mask2D.write_image('mask2D_test.hed')
	return smask2D

#==========================
def createMask(params):
	nx = params['nx']
	apix = params['apix']

	img = EMData(nx,nx)
	img.to_zero()

	#add 3 degrees to overlap with the neighboring density
	overlap=3*math.pi/180.0
	alpha = math.pi/2 - math.pi/params['pf'] - overlap
	#creating wedge-slice
	for x,y in ((x,y) for x in range(0,nx) for y in range(nx/2,nx)):
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
		
	smask2D = circularMask2D(nx)
	img.mult(smask2D)
	wedge = EMData(nx,nx,nx)
	twist = params['twist']
	rise = params['rise']
	alpha = 360+(params['pf']*twist)

	for z in range(nx):
		l = params['pf']*rise
                rot = alpha/l*apix
                finalrot = ((z-nx/2)*rot)/3
                finalrot = locateGoodPF(params)
                if params['pf'] == 14:
                        finalrot = finalrot + 3
                t = Transform()
                t.set_rotation({"type":"2d","alpha":-finalrot})
                newslice = img.process("xform",{"transform":t})
                wedge.insert_clip(newslice,(0,0,z))

        if params["decor"] == "kinesin":
                print "decor = kinesin"
                ymsk = int(148/apix)
                xmsk = int(49/apix)
                mskrad = int(20/apix)
                if params['pf'] == 12:
                        ymsk = int(142/apix)
                        xmsk = int(42/apix)
                        mskrad = int(16/apix)
                        
                # see if mask is near the edge:
                edge=ymsk*math.atan(math.pi/params['pf'])
                if (abs(xmsk)+mskrad)>=edge:
                        # distance for corresponding positive mask
                        edge = int(2*edge)
                        xmsk2 = int(math.copysign(edge-abs(xmsk),xmsk)*-1)
                        # take max of 1 mask
                        avgr = Averagers.get("minmax",{"max":1})
                        avgr.add_image_list([wedge,kinesinMask(nx,mskrad,xmsk2,ymsk,pos=True)])
                        wedge=avgr.finish()
                # multiply 0 mask
                wedge *= kinesinMask(nx,mskrad,xmsk,ymsk)

	if params['opposite'] is True:
		t = Transform({"type":"spider","psi":180})
		wedge.process_inplace("xform",{"transform":t})

	return wedge

#===========================
def applySeamSym(params):
	"""
	apply seam symmetry based on results from Egelman search
	"""
	vol = params['vol']
	rise = params['rise']
	rot = params['twist']
	apix = params['apix']

	start_time = time.time()

	# find protofilament number from rotation
	sym=int(round(360.0/abs(rot)))

	rise/=apix
	# apply protofilament symmetry
	sumvol = vol.copy()
	pfoffset=int(sym/2)
	for pnum in range(-pfoffset,sym-pfoffset):
		if pnum==0: 
			continue
		if params["decor"]=="EB" and pnum == -pfoffset:
			continue
		ang = rot*pnum
		trans = -(rise*pnum)
		t = Transform({"type":"spider","psi":ang})
		t.set_trans(0,0,trans)
		volcopy = vol.process("xform",{"transform":t})
		sumvol.add(volcopy)

	sumvol.process_inplace("normalize")
	sumvol.write_image('volOverSym.mrc')
	print "Real-space symmetry applied in %.2f minutes"%((time.time()-start_time)/60.0)

	return sumvol

#===========================
def regenerateFromPF(params,wedgemask):
	"""
	mask out one away from seam, regenerate microtubule with seam 
	"""
	import shutil,subprocess

	# convert rise to pixels
	nx = params['nx'] 
	rise = params['rise']/params['apix']
	twist = params['twist']

	# apply protofilament symmetry
	#sumvol = params['vol']*wedgemask
	
	sumvol = EMData(nx,nx,nx)
	sumvol.to_zero()
	#sumvol.write_image("pf.mrc")
	pfoffset=int(params['pf']/2)
	if params['opposite'] is True:
		pfoffset=0

	start_time = time.time()

	for pnum in range(-pfoffset,params['pf']-pfoffset):
		#if pnum == 0:
		#	continue
		print "preparing copy %i"%pnum
		ang = twist*pnum
		trans = -(rise*pnum)
		ang*=-1
		trans*=-1
		t = Transform({"type":"spider","psi":ang})
		t.set_trans(0,0,trans)
		volcopy = params['vol'].process("xform",{"transform":t})
		seammaskcopy = wedgemask.process("xform",{"transform":t})
		#seammaskcopy.write_image("wedgemask_%d.mrc"%pnum)
		#svol=sumvol*(1-seammaskcopy)
		#vcopy=volcopy*seammaskcopy

		sumvol = sumvol*(1-seammaskcopy)+volcopy*seammaskcopy

	#sumvol.process_inplace("xform",{"transform":t})

	print "Seamed MT regenerated in %.2f minutes"%((time.time()-start_time)/60.0)
	if params['vrise'] is not None and params['vtwist'] is not None:
		start_time = time.time()
		edgeMask(params)
		sumvol1 = sumvol.helicise(params['apix'],params['vrise'],params['vtwist'])
		# because symmetry is only applied in one direction, flip & apply again
		t = Transform({"type":"spider","theta":180})
		sumvol.process_inplace("xform",{"transform":t})
		sumvol2 = sumvol.helicise(params['apix'],params['vrise'],params['vtwist'])
		sumvol2.process_inplace("xform",{"transform":t})
		sumvol = sumvol1+sumvol2
		del sumvol1,sumvol2
		print "Dimer symmetry applied in %.2f minutes"%((time.time()-start_time)/60.0)

	params['vol']=sumvol.process("normalize")
	

#===========================
def kinesinMask(nx,rad,cx,cy,pos=False):
	# soft edge cylinder mask for kinesin position
	img = EMData(nx,nx)
	img.to_one()
	if pos is True:
		img.to_zero()

	# outer radius
	orad = (rad+rad*.5)

	if abs(cy) > (nx/2-orad) : cy = int((cy/abs(cy))*(nx/2-orad))
	if abs(cx) > (nx/2-orad) : cx = int((cx/abs(cx))*(nx/2-orad))
	for x,y in ((x,y) for x in range(-nx/2,nx/2) for y in range(-nx/2,nx/2)):
		r2 = x**2+y**2
		if r2 < (orad*orad):
			if r2 < rad*rad:
				val=1
			else:
				diff=orad**2-rad**2
				val=1-((r2-rad*rad)/(diff))
			if pos is True:
				img.set(nx/2-x+cx,nx/2+y+cy,val)
			else:
				img.set(nx/2+x+cx,nx/2+y+cy,1-val)
		
	cylmask = EMData(nx,nx,nx)
	twist = params['twist']
	rise = params['rise']
	alpha = 360+(params['pf']*twist)
	for z in range(nx):
		l = params['pf']*rise
		rot = alpha/l*params['apix']
		finalrot = (z-nx/2)*rot
                #finalrot = params['finalrot']
                finalrot = locateGoodPF(params)
                if params['pf'] == 14:
                        finalrot = finalrot + 3
		t=Transform()
		t.set_rotation({"type":"2d","alpha":-finalrot}) #commented out finalrot/3
		newslice=img.process("xform",{"transform":t})
		cylmask.insert_clip(newslice,(0,0,z))
        cylmask.write_image("kinesinmask.mrc")
	return cylmask

#==========================
def edgeMask(params):
	"""
	create a 3D cylinder mask to remove edges and artifacts from symmetrization
	"""
	nx = params['nx']
	nxm = int(nx+(nx*0.3))
	nym = int(nx-(nx*0.3))

	params['vol'].write_image("tmp.mrc")

	mask2d = EMData(nxm,nym)
	mask2d.to_one()
	mask2d.process_inplace("mask.decayedge2d",{"width":nx*0.1})
	mask2d.clip_inplace(Region(int(nx*0.3)/2,-int(nx*0.3)/2,nx,nx))
	mask=EMData(nx,nx,nx)
	for i in xrange(nx):
		mask.insert_clip(mask2d,(0,0,i))
	
	t = Transform({"type":"spider","theta":90.0,"phi":90.0})
	mask.process_inplace("xform",{"transform":t})

	irad = int(nx/2*0.85)
	orad = (nx/2)-2
	falloff = orad - irad
	for i in xrange(nx):
		slice2d = mask.get_clip(Region(0,0,i,nx,nx,1))
		slice2d.process_inplace("mask.gaussian",{"inner_radius":irad,"outer_radius":falloff})
		mask.insert_clip(slice2d,[0,0,i])
	params['vol']*=mask

#===========================
def CreateEBSeamMask(params):
        nx = params['nx']
        rise = params['rise']
        twist = params['twist']
        apix = params['apix']
        # soft-edged cylinder mask to remove EB at the seam
        sym=int(round(360.0/abs(twist)))
        alpha = 360+(sym*twist)
        lrise = sym*rise
        rot = alpha/lrise*apix
        # hard coded for the moment
        #seammask = kinesinMask(nx,5,10,-50,rot)
        seammask = kinesinMask(nx,int(14/apix),int(28/apix),int(-137/apix))
        return seammask
       
       
# Generate 2D slices to be inserted into mask3D volume
def createMask2D(params):
	from itertools import product
	apix = params['apix']
	orad = float(params['orad'])/apix
	irad = float(params['irad'])/apix
	nx = params['nx']
	falloff_r = 30			# use steeper falloff
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
	if params['zrad']:
		zrad = float(params['zrad'])/apix
	else:
		zrad = int(nx/2*0.8)
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
if __name__ == "__main__":
	#getEMANPath()
	from EMAN2 import *
	from sparx import *
	params=setupParserOptions()
	checkConflicts(params)
	# normalize to mask
#	params['vol']+=0.03407
#	params['vol'].process_inplace("filter.flattenbackground",{"radius":params['nx']/2-4})
#	params['vol'].process_inplace("normalize.circlemean")
#	bg = params['vol'].get_circle_mean()
#	print "background:",bg
#	bg = params['vol'].get_edge_mean()
#	print "background:",bg
#	params['vol']-=params['vol'].get_circle_mean()
	#params['vol'].write_image("tmp.mrc")
	foundWedgeMask = False
	if params['applysym'] is True:
		params['vol'] = applySeamSym(params)
		params['vol'].write_image('output_applySym.mrc')
	
	try: 
		wedgemask = EMData('wedgemask.mrc')
		params['finalrot']=-1
		foundWedgeMask = True
	except:
		print "cannot find existing wedgemask.mrc, will create a new one"
	if not foundWedgeMask:
		wedgemask = createMask(params)
		wedgemask.write_image("wedgemask.mrc")

	regenerateFromPF(params,wedgemask)
	
	if params["decor"] == "EB":
		seammask = CreateEBSeamMask(params)
		params['vol']*= seammask
	
	mask2D = createMask2D(params)
	mask3D = createMask3D(params,mask2D)
	mask3D.write_image('mask3D.mrc')
	params['vol']*= mask3D
	
#	edgeMask(params)
	params['vol'].write_image("output.mrc")
