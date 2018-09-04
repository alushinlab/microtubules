#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will calculate the widths of MT class averages based on the temp.spi outputs from InspectMTStack_norm_measurePeaks.spi
#
#import modules, this could probably stand to be pruned
#

# import system / house keeping modules
import os
import subprocess
import sys
import optparse
import numpy
import scipy.signal
#
#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -i <file_list> --min <xmin> --max <xmax> --apix <angstroms per pixel> --val1 <val1> --val2 <val2>")
	parser.add_option("-i",dest="input",type="string",metavar="FILE",
		help="list of temp files to quantify")
	parser.add_option("--min",dest="min",type="int",metavar="INT",
		help="minimum x value for peak search in pixels")
	parser.add_option("--max",dest="max",type="int",metavar="INT",
		help="maximum x value for peak search in pixels")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
		help="Angstroms per pixel for width calculation")
	parser.add_option("--val1",dest="val1",type="int",metavar="int",default=50,
		help="Peak width 1.  Should not adjust unless there is an error.  Default 50")
	parser.add_option("--val2",dest="val2",type="int",metavar="int",default=50,
		help="Peak width 2.  Should not adjust unless there is an error.  Default 50")
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 5:
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
#=============================
	
if __name__ == "__main__":
	params=setupParserOptions()
	infile=params['input']
	checkExistence(infile)
	
	#set up variables
	minimum=params['min']
	maximum=params['max']
	apix=params['apix']
	val1=params['val1']
	val2=params['val2']

	#see if program has run already
	if os.path.exists('distances.txt'):
		print "Warning, this program has already been run in this directory.  Remove previous outputs and run again."
		sys.exit()

	#make a tuple of fitting widths
	widths=val1,val2

	#let's quantify microtubules!
	with open(infile,'r') as a:
		file_list=a.readlines()
	for filename in file_list:
		with open (filename.strip(),'r') as a:
			values_withjunk = a.readlines()
		values_withjunk.pop(0) # pop first line from spider file, is a comment
		values = []
		for value in values_withjunk:
			value_split=value.split()
			values.append(value_split[2])	#get third column
		values_array = numpy.array(values,dtype=float) #convert values from string to an array of floats
		relevant_values = values_array[minimum:maximum] #only consider values within the designated range
		peaks=scipy.signal.find_peaks_cwt(relevant_values,widths) #find peaks
		print peaks
		peaks=peaks+minimum #Scale by pixel size and account for the offset due to boundaries.  Numpy arrays are sweet!
		peaks_scaled=apix*peaks
		if len(peaks_scaled) != 2:
			print "Warning: more or less than 2 peaks found for file %s.  The distance calculation for this file is not reliable"%(filename.strip())
			peaks=['NaN','NaN']
			peaks_scaled=['NaN','NaN']
			distance='NaN'
		else:
			distance=abs(peaks_scaled[1]-peaks_scaled[0])
		with open ('distances.txt','a') as distance_out:
			distance_out.write("%s\n"%str(distance))
		with open ('peaks.txt','a') as peaks_out:
			peaks_out.write("%s     %s     %s     %s\n"%(str(peaks[0]),str(peaks[1]),str(peaks_scaled[0]), str(peaks_scaled[1])))
					
	
			
	

	
