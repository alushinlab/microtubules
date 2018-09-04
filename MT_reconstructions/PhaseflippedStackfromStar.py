#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will combine particles from RELION substacks into a single phase-flipped stack 
# suitable for EMAN2 / SPARX IHRSR In the order indicated in a STAR file
#
#
# 
# import system / house keeping modules
import os
import subprocess
import sys
import optparse
# import EMAN2 and sparx python bindings for image handling
# EMAN2 is required to run this script
from EMAN2 import *
from sparx import *
#
#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -i <input STAR> -o <output stack>")
	parser.add_option("-i",dest="star",type="string",metavar="FILE",
		help="input star file")
	parser.add_option("-o",dest="outstack",type="string",metavar="FILE",default='start.img',
		help="ouput stack in imagic, default = start.img")
	
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
def checkExistence(f_name):
	if not os.path.exists(f_name):
		print "Error: file %s does not exist.  Exiting..."%f_name
		sys.exit()
#=========================
# Adapted from John Rubinstein's starfile handler
def find_index(label,labels):
	for i, s in enumerate (labels):
		if label in s:
			return i
	return False
#======================
def Split_header(name):
	#find the last line that starts with _, this is the end of the header
	#also find star file column labels
	f=open(name,'r')	
	lines=f.readlines()
	f.close()
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
def get_imageinfo(star_line,index):
	splitline=star_line.split()
	#split at @ to get particle position and file
	image_record = splitline[index]
	record=image_record.split('@')
	relion_index,path=record[0],record[1]
	e2_index=int(relion_index)-1
	
	return e2_index,path
	
#======================	
def StripBlankLinesFromEnd(mylist):
	while mylist[-1] == ' \n':
		mylist.pop(-1)
	return mylist
#=============================


if __name__=='__main__':
	params=setupParserOptions()    
	star=params['star']	
	outstack=params['outstack']
	
	checkExistence(star)

	#split up our star file into header, body, and labels
	starheader,starbody,starlabels=Split_header(star)

	#remove irritating trailing newlines from star body:
	starbody=StripBlankLinesFromEnd(starbody)
	
	#Find the relevant indexes in the star file
	imagecolumn = find_index('_rlnImageName',starlabels)
	
	#CTF	
	df1	= find_index('_rlnDefocusU',starlabels)
	df2	= find_index('_rlnDefocusV',starlabels)
	kV	= find_index('_rlnVoltage',starlabels)
	Cs	= find_index('_rlnSphericalAberration',starlabels)
	Apix	= find_index('_rlnDetectorPixelSize',starlabels)
	AmpC	= find_index('_rlnAmplitudeContrast',starlabels)


	#For simplicity we will ignore astigmatism.  This is what we always did
	#In appion

	AngAst	 = 0.0 
	AstigAmp = 0.0
		

	#initate counter and EMData() objects for input and output stacks
	Stack_index=0	
	particle=EMData()
	flipped_particle=EMData()
	
	#Now iterate through body to build our new stack
	for line in starbody:
		in_index,in_path=get_imageinfo(line,imagecolumn)
		particle.read_image(in_path,in_index)
		
		entries	= line.split()		
		defocus = (float(entries[df1]) + float(entries[df2]))/2
		AmpCont = 100*float(entries[AmpC])
		
		#generate CTF object and phase flip particle		
		ctf = generate_ctf([defocus,entries[Cs],entries[kV],entries[Apix],0,AmpCont,AstigAmp,AngAst])
		flipped_particle=filt_ctf(particle,ctf,True,1,1) 

		
		flipped_particle.write_image(outstack,Stack_index)
		Stack_index += 1
	
	
	
