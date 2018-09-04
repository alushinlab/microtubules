#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu
# Modified by Matthew Reynolds
# This script is meant to be a companion to a shell script that parallelizes the 
# Gain_normalize_and_bin_movies task. It requires input of the c-shell file and number
# of partitions to use for parallelization. 
#
# Requires: .csh file in the same directory, a defects.txt file, and a reference.mrc file
# Input Requirements: a valid .csh file must be entered in the --shProg input, and a non-
# negative integer must be entered for the --numGroups input
#
# import modules
#
import os
import subprocess
import sys
import optparse
import time

#
#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog --shProg <shell Script to be parallelized> --numGroups <number of partitions to make during parallelization>")
	parser.add_option("--shProg",dest="shProg",type="string",metavar="string",
		help="Input the shell program that you want to parallelize")
	parser.add_option("--numGroups",dest="numGroups",type="int",metavar="int",
		help="Number of partitions to make during parallelization")

	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	# if improper number of inputs entered
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params
#=========================
# checks the existance of the .csh program in the current directory
def checkExistence(f_name):
	if not os.path.exists(f_name):
		print "Error: file %s does not exist.  Exiting..."%f_name
		sys.exit()
#=========================
# checks that the group number entered is a positive integer for proper parallelization
def ensurePositiveInt(groupNum):
	if groupNum < 1:
		print "Error: group number %s is a non-positive integer.  Exiting..."%groupNum
		sys.exit()
#===========================
if __name__ == "__main__":
	params=setupParserOptions()
	# Make sure the shell script exists in the same directory and numGroups > 0
	checkExistence(params['shProg'])
	ensurePositiveInt(params['numGroups'])
	
	# determines number of files in the current working directory that end with .tif
	files=os.listdir(os.getcwd())
	files_tif = [i for i in files if i.endswith('.tif')]
	numFiles = len(files_tif)

	# determines the cutoff points for each file grouping
	partitions = [((numFiles//params['numGroups']) * i) for i in range(params['numGroups'])]

	# if only one group, then run in series	
	if params['numGroups'] == 1:
		prepCMD1="./%s %s %s > %s.log"%(str(params['shProg']),str(0),str(numFiles),str(params['shProg'][:-4]+str(0)))
		subprocess.Popen(prepCMD1,shell=True)
	
	# if more than one group, then run in parallel
	else:
		k=0
		for j in partitions[:-1]:
			prepCMD1="./%s %s %s > %s.log"%(str(params['shProg']),str(j),str(j+partitions[1]),str(params['shProg'][:-4]+str(k)))
			subprocess.Popen(prepCMD1,shell=True)
			k=k+1
			time.sleep(1)
		# treat last element differently because the partitions may not have exactly even entries
		# extras are put into the last group
		lastElement = partitions[-1]
		prepCMD1="./%s %s %s > %s.log"%(str(params['shProg']),str(lastElement),str(lastElement+partitions[1]*2),str(params['shProg'][:-4]+str(k)))
		subprocess.Popen(prepCMD1,shell=True)
		time.sleep(1)


