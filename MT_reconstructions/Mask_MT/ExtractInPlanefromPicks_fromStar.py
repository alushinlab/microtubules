#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will extract relion PSI tilt prior from picks, then output a spider doc file for applying a softmask
# This is useful for soft-masking helical segments prior to running EMAN2 reference-based classification / reconstruction
#
# 
# import system / house keeping modules
import os
import subprocess
import sys
import optparse

#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -i <input STAR>")
	parser.add_option("-i",dest="star",type="string",metavar="FILE",
		help="input star file")
	
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 2:
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
#======================
def get_record(star_line,index):
	splitline=star_line.split()
	record=splitline[index]
	return record	
#======================	
def StripBlankLinesFromEnd(mylist):
	while mylist[-1] == ' \n':
		mylist.pop(-1)
	return mylist
#=============================


if __name__=='__main__':
	params=setupParserOptions()    
	star=params['star']	
	checkExistence(star)
	if os.path.exists("inplane_from_star.spi"):
		print "Warning: this script has alread been run in this directory.  Delete inplane_from_star.spi and run again if desired."
		sys.exit()
	
	#split up our star file into header, body, and labels
	starheader,starbody,starlabels=Split_header(star)

	#remove irritating trailing newlines from star body:
	starbody=StripBlankLinesFromEnd(starbody)
	
	#Find the index of Psi prior from picks in the star file
	imagecolumn = find_index('_rlnAnglePsiPrior',starlabels)

	#initate counter, spider index begins with 1
	Stack_index=1	
	#Now iterate through body to build our new stack
	for line in starbody:
		angle=get_record(line,imagecolumn)
		angle=(float(angle))+90		
		outline="%5d 1%11.4f\n"%(int(Stack_index),float(angle))
		with open("inplane_from_star.spi","a") as f:
			f.write(outline)	 		
		Stack_index += 1
	

	
