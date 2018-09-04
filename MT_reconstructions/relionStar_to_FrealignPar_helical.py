#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will convert a relion STAR parameter file into a Frealign .par file.
# For simplicity it will produce a file in Frealign 8 format, which is still recognized by Frealign 9  
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
	parser.set_usage("%prog -i <input STAR> -o <output .par>")
	parser.add_option("-i",dest="star",type="string",metavar="FILE",
		help="input star file")
	parser.add_option("-o",dest="outpar",type="string",metavar="FILE",
		help="ouput .par FREALIGN parameter file")
	
	
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
def Par_Line(count,psi,theta,phi,shx,shy,mag,film,df1,df2,angast):
	#String formatting for par file entry taken from EMAN2FREALIGN.py, Nogales lab
	parline=("%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.1f%6d%9.1f%9.1f%8.2f%7.2f%8.2f\n" %(count,psi,theta,phi,shx,shy,mag,film,df1,df2,angast,0,0))
	return parline
#======================
def StripBlankLinesFromEnd(mylist):
	while mylist[-1] == ' \n':
		mylist.pop(-1)
	return mylist
#=============================

if __name__=='__main__':
	params=setupParserOptions()    
	star=params['star']	
	outpar=params['outpar']
	
	checkExistence(star)
	
	#split up our star file into header, body, and labels
	starheader,starbody,starlabels=Split_header(star)
	
	#remove irritating trailing newlines from star body:
	starbody=StripBlankLinesFromEnd(starbody)

	#Find the indexes of all the parameters we need for our FREALIGN output file
	psi	= find_index('_rlnAnglePsi',starlabels)
	theta	= find_index('_rlnAngleTilt',starlabels)
	phi	= find_index('_rlnAngleRot',starlabels)
	shx	= find_index('_rlnOriginX',starlabels)
	shy	= find_index('_rlnOriginY',starlabels)			
	df1	= find_index('_rlnDefocusU',starlabels)
	df2	= find_index('_rlnDefocusV',starlabels)
	angast	= find_index('_rlnDefocusAngle',starlabels)
	
	#Relion helical tube ID is counted by micrograph.  So we need to calculate it for FREALIGN
	helID	= find_index('_rlnHelicalTubeID',starlabels)
	#Always just set mag to equal 10000, per convention
	mag	= 10000
	
	#Now loop over body entries to make our Frealign par file
	Count  	= 1
	Filament_count = 1
	Current_fil    = 1
	
	with open(outpar,'w') as f:
		for line in starbody:
			entries=line.split()
			new_fil=int(entries[helID])
			if new_fil != Current_fil:
				Filament_count += 1
				Current_fil = new_fil
			ParLine=Par_Line(Count,float(entries[psi]),float(entries[theta]),float(entries[phi]),-1*float(entries[shx]),-1*float(entries[shy]),
				mag,Filament_count,float(entries[df1]),float(entries[df2]),float(entries[angast]))
			f.write(ParLine)
			Count += 1

	
