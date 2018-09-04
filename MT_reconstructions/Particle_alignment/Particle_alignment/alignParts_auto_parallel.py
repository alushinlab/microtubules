#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will prepare inputs for a parallel run of alignparts_lmbfgs, as well as a star file pointing to outputs
# Then launch N parallel processes
# Be careful with N: stay within memory constraints and the read-write capabilities of disk

#
#import modules, this could probably stand to be pruned
#

# import system / house keeping modules
import os
import subprocess
import sys
import optparse

#
#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -i <input star file> -m <movies directory> --box <box size> -t <alignparts_template> -N <number of processes to launch>")
	parser.add_option("-i",dest="instar",type="string",metavar="FILE",
		help="input particles star file")
	parser.add_option("-m",dest="movie_dir",type="string",metavar="string",
		help="directory containing gain-normalized movies")
	parser.add_option("--nx",dest="nx",type="int",metavar="int",default=3838,
		help="Size of movie frames in x, default 3838")
	parser.add_option("--ny",dest="ny",type="int",metavar="int",default=3710,
		help="Size of movie frames in x, default 3710")
	parser.add_option("--box",dest="box_size",type="int",metavar="int",
		help="Box size of particles")
	parser.add_option("-t",dest="template",type="string",metavar="FILE",default="alignparts_lmbfgs_template.csh",
		help="template c-shell script to run alignparts.  Must edit to set right parameters for your dataset")
	parser.add_option("-N",dest="Nprocs",type="int",metavar="int",
		help="Number of instances of align_parts to launch.  Be cautious.")	
	
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 8:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params
#=========================
def checkExistence(f_name):
	if not os.path.exists(f_name):
		print "Error: path %s does not exist.  Exiting..."%f_name
		sys.exit()
#=========================
#taken directly from stackoverflow user Max Shawabkeh
#http://stackoverflow.com/questions/2130016/splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))
#========================
def modify_string(old,find,replace):
	new=[]
	for line in old:	
		new.append(line.replace(find,replace))
	return new
#======================
def Split_header(name):
	#find the last line that starts with _, this is the end of the header
	#also find star file column labels
	with open(name,'r') as f:	
		lines=f.readlines()
	labels=[]
	for line in lines:
		if line.startswith('_'):
			labels.append(line)
			endheader=lines.index(line)
	#now split the list into 2
	endheader+=1
	print "Header ends at line %s"%str(endheader)	
	header=lines[0:endheader]
	body=lines[endheader:]
	return header,body,labels

#======================
def StripBlankLinesFromEnd(mylist):
	while mylist[-1] == ' \n':
		mylist.pop(-1)
	return mylist
#=========================
# Adapted from John Rubinstein's starfile handler
def find_index(label,labels):
	for i, s in enumerate (labels):
		if label in s:
			return i
	return False
#
def checkParticle(partX,partY,imX,imY,boxsize):
	xposs=partX-boxsize/2
	xposf=partX+boxsize/2-1
	yposs=partY-boxsize/2
	yposf=partY+boxsize/2-1
	if xposs<0 or xposf>imX or yposs<0 or yposf>imY:
		return False
	return True
#taken directly from starfile handler:
def writestarline(outfile,records):
    # JLR 7/16
    """Write a record (line) to an already open starfile"""
    for item in records:
        outfile.write(item+'  ')
    outfile.write('\n')

#========================
#functions for reordering star file body if necessary
#
#little function to make a key pair from the appropriate labels
def makekey(a,index1,index2):
	val1=a.split()[index1]
	val2=a.split()[index2]
	pair=(val1,val2)
	return pair
# reorder function
def reorderStar(inlist,Labels):
	micro=find_index('_rlnMicrographName',Labels)
	print "Micrograph column is %s"%str(micro)
	image=find_index('_rlnImageName',Labels)
	print "Image column is %s"%str(image)

	#now use nifty little function to sort the body
	outlist=sorted(inlist,key = lambda x : makekey(x,micro,image))
	
	return outlist


#=========================
if __name__ == "__main__":
	params=setupParserOptions()

	#gather inputs and verify their existence
	instar=params['instar']
	template=params['template']
	Nprocs=params['Nprocs']
	movie_dir=params['movie_dir']
	nx=params['nx']
	ny=params['ny']
	box=params['box_size']

	checkExistence(instar)
	checkExistence(template)
	checkExistence(movie_dir)
	
	if box is None:
		print "Please indicate box size with --box flag"
		sys.exit()

	if Nprocs is None:
		print "Please indicate number of processes with -N flag"
		sys.exit()
		
	#Split input star into components:
	header,body,labels=Split_header(instar)
	body=StripBlankLinesFromEnd(body)

	#Check if star file needs reordering for alignparts to work, and warn user if so:
	body_reordered=reorderStar(body,labels)
	if body_reordered != body:
		print "Warning: input star file was out of order. It has been reordered to ensure proper function of alignparts_lmbfgs"
		body = body_reordered
	
	
	#Find indexes of relevant metadata
	coordX_index	= find_index('_rlnCoordinateX',labels)
	coordY_index	= find_index('_rlnCoordinateY',labels)
	micro_index	= find_index('_rlnMicrographName',labels)
	image_index	= find_index('_rlnImageName',labels)
	
	outstar='%s_lmbfgs.star'%(os.path.splitext(instar)[0])
	movie_list=[]
	
	##################################################################	
	#This section adapted from John Rubinstein's starfile handler
	
	#intialize counters
	totparts = 0; nmovies = 0; nskipped = 0; nparts = 0 ; prevmic = 'x'; firstmic = True; nskipped_micro=0

	#Now iterate through body to generate coordinate files and output star file
	with open(outstar,'w') as Outstar:
		for line in header:
			Outstar.write(line)		
		for line in body:
			#unpack record and advance counters.
			record=line.split()
			totparts+=1
                	nparts+=1 # particle number for this movie (even if excluded)
			
			#Get particle coordinates and check if it is in bounds			
			coordX=int(float(record[coordX_index]))
			coordY=int(float(record[coordY_index]))
			partStatus = checkParticle(coordX,coordY,nx,ny,box)

			if partStatus is False:
				print 'SKIPPED particle {:05} in movie {:04}: Coordinate out-of-bounds for specified boxsize'.format(nparts,nmovies)
                   		nskipped +=1
		    		nskipped_micro +=1						
			else:
                    		if record[micro_index] != prevmic:
                        		if not firstmic:
                           			outcoordfile.close()
                      			else:
                            			print "First micrograph"
                            			firstmic = False
                       			prevmic = record[micro_index]
				
					#write out coordfile name in movies directory.  Setting this up is a bit tedious
					outcoordbasename=os.path.splitext(os.path.basename(record[micro_index]))[0]	
					outcoordpath=os.path.join(movie_dir,'%s.coord'%(outcoordbasename))

					#append parter movie to movie list for easy adaptation of alignparts_parallel					
					outcoord_moviepartner='%s.mrcs\n'%(outcoordbasename)		
					movie_list.append(outcoord_moviepartner)

                        		outcoordfile = open(outcoordpath,'w')
                        		nparts=0
					nskipped_micro=0
                        		nmovies+=1
                        		outcoordfile.write('        x         y      density\n')
                    		coordline = "{0:>10} {1:>9} {2:>12}".format(coordX,coordY,'0.0\n')
		   		outcoordfile.write(coordline) # write line to the coordfile
		    		#correct image name record for output star file
		    		npartsincluded=nparts-nskipped_micro+1
		    		stack_filename=record[image_index].split('/')[-1]
		    		record[image_index]='%s@%s/%s'%(str(npartsincluded).zfill(6),'Particles_lmbfgs',stack_filename)		  
                    		writestarline(Outstar,record) # write line to the starfile
        	print "Total number of particles recognized",totparts
       		print "Number of skipped particle images",nskipped
	#############################################################################
	### This section will now prepare inputs for a parallel run of alignparts and run it
	

	if os.path.exists('Particles_lmbfgs'):
		print "Error, Particles_lmbfgs directory already exists, perhaps you have previously run this script.  Please move or delete it, and restart."
		sys.exit()

	os.mkdir('Particles_lmbfgs')
		
	with open (template,'r') as f:
		Template=f.readlines()

	#split movies list into Nprocs ~equal parts:
	movieLists=list(split(movie_list,Nprocs))

	#loop through split list to prepare input files for each parallel run of alignparts
	for i in xrange(Nprocs):
		#write a movie list file for process i		
		currentmlist='movies_%s.txt'%(str(i))			
		with open(currentmlist,'w') as f:		
			for line in movieLists[i]:
				f.write(line)

		#make a list of coord files with our nifty modify string module
		clist=modify_string(movieLists[i],'.mrcs','.coord')		
		currentclist='coords_%s.txt'%(str(i))
		with open(currentclist,'w') as f:
			for line in clist:
				f.write(line)
		
		#modify and write out appropriate alignparts_script
		currentscript='alignparts_script_%s.csh'%(str(i))
		script_rightmovies=modify_string(Template,'movies.txt',currentmlist)
		script_rightcoords=modify_string(script_rightmovies,'coords.txt',currentclist)
		
		with open (currentscript,'w') as f:
			for line in script_rightcoords:
				f.write(line)
		
		#make script executable
		cmd = 'chmod 755 %s'%currentscript
		subprocess.Popen(cmd,shell=True).wait()	
	#now execute Nprocs scripts
	for i in xrange(Nprocs):
		currentscript='alignparts_script_%s.csh'%(str(i))
		cmd ="./%s > %s.log"%(currentscript,currentscript[:-4])
		print cmd
		subprocess.Popen(cmd,shell=True)
	
		
