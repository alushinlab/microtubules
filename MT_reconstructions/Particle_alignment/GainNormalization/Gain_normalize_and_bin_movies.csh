#!/bin/csh
# Script to generate doseweighted and nondoseweighted drift-corrected images with unblur from compressed tifs
# Requires tifs, defects file, and a gain reference, appropriately re-oriented (flipped, rotated) in the same directory
#
# Edit the following variables to set up parameters appropriate for your experiment
#
#
####################################################

set defects=defects_deltaTMD0TAPpeptide_0001.txt
set reference=reference_fliprot90.mrc
set apix=0.5
#Note that apix should be the superres pixel size
set num_frames=50
set dose_rate=1.6
set kV=300


set dir1=Movies_Bin2
@ i = $1
@ j = $2
@ k = 0

######################################################

# Make directory to store movies
if (! -d $dir1) then
	mkdir $dir1
endif

# Now gain normalize and bin each file between the specified values
foreach file (*.tif)
if ( $k >= $i ) then
	if ( $k < $j) then
		echo $i $j $k
		echo $file
		#do operations:
		#first gain normalize the movie
		clip mult -m 2 -D $defects $file $reference ${file:r}_normalized.mrc
		#bin by 2 since we collect superres movies
# Watch out, this is whitespace sensitive
resample.exe <<eot
${file:r}_normalized.mrc
${file:r}_aligned_bin2.mrc
NO
NO
3710
3838
eot

		#move file to save, delete bin 0.5 file
		mv ${file:r}_aligned_bin2.mrc Movies_Bin2/${file:r}_aligned_bin2.mrcs
		rm -f ${file:r}_normalized.mrc
	endif
endif
@ k += 1

end
