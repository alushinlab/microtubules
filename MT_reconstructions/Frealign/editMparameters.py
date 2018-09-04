#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will edit Frealign mparameters files to update symmetry and refinment round number
# It is meant to interface with refineMTfre_v9p11 scripts and mparameters_template file

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
	parser.set_usage("%prog -i <input mparameter file> -o <output mparameter file> --asym <symmetry setting> --num <round of refinement> --msk <parameter mask> --aS <asym string> --fS <first round number string>  --lS <last round number string> --mS <parameter mask string>")
	parser.add_option("-i",dest="in_Mpar",type="string",metavar="FILE",default="mparameters_template",
		help="input mparameter file, default mparameters_template")
	parser.add_option("-o",dest="out_Mpar",type="string",metavar="FILE",default="mparameters",
		help="output mparameter file, default mparameters")
	parser.add_option("--asym",dest="asym",type="string",metavar="string",
		help="Symmetry setting to update, sane options are 0 (corresponding to C1) or HP")
	parser.add_option("--num",dest="num",type="int",metavar="int",
		help="Round of refinement")
	parser.add_option("--msk",dest="msk",type="int",metavar="int",default="1",
		help="Parameter mask value.  0 disables shift/angle modification, 1 enables.  0 is useful for priming frealign.  default 1")
	parser.add_option("--aS",dest="aS",type="string",metavar="string",default="Symmetry              X",
		help="Symmetery string to replace, default Symmetry              X")
	parser.add_option("--fS",dest="fS",type="string",metavar="string",default="start_process         X",
		help="Symmetery string to replace, default start_process         X")
	parser.add_option("--lS",dest="lS",type="string",metavar="string",default="end_process           X",
		help="Symmetery string to replace, default end_process           X")	
	parser.add_option("--mS",dest="mS",type="string",metavar="string",default="X X X X X",
		help="Mask string to replace, default X X X X X")
	
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
		print "Error: path %s does not exist.  Exiting..."%f_name
		sys.exit()
#========================
def modify_string(old,find,replace):
	new=[]
	for line in old:	
		new.append(line.replace(find,replace))
	return new
#=========================
if __name__ == "__main__":
	params=setupParserOptions()
	
	#check for existence of mparameters template
	checkExistence(params['in_Mpar'])

	#set up strings to make edits
	new_Sym="Symmetry              %s"%(params['asym'])
	new_First="start_process         %s"%(str(params['num']))
	new_Last="end_process           %s"%(str(params['num']))
	new_Mask="%s %s %s %s %s"%(str(params['msk']),str(params['msk']),str(params['msk']),str(params['msk']),str(params['msk']))

	#find and replace
	with open(params['in_Mpar'],'r') as f:
		in_Mpar = f.readlines()
	
	temp=modify_string(in_Mpar,params['aS'],new_Sym)
	temp=modify_string(temp,params['fS'],new_First)
	temp=modify_string(temp,params['lS'],new_Last)
	temp=modify_string(temp,params['mS'],new_Mask)
	
	with open(params['out_Mpar'],'w') as f:
		for line in temp:
			f.write(line)


	
	

