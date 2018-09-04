#!/bin/csh -f
#
#   Control script to submit multiple jobs on a cluster using
#   the Sun Grid Engine. Each job processes N particles. N is
#   specified in the 'mparameters' file as 'increment'.
#
cp mparameters_halves_helical mparameters_run
set stack_a = `grep stack_a mparameters_run | awk '{print $2}'`
set stack_b = `grep stack_b mparameters_run | awk '{print $2}'`
set stack_m = `grep stack_m mparameters_run | awk '{print $2}'`

set model_a = `grep model_a mparameters_run | awk '{print $2}'`
set model_b = `grep model_b mparameters_run | awk '{print $2}'`
set model_m = `grep model_m mparameters_run | awk '{print $2}'`

set start = `grep start_process mparameters_run | awk '{print $2}'`
set end   = `grep end_process mparameters_run | awk '{print $2}'`
set first = `grep first_particle mparameters_run | awk '{print $2}'`
set last  = `grep last_particle mparameters_run | awk '{print $2}'`
set incr  = `grep increment mparameters_run | awk '{print $2}'`
set mode = `grep MODE mparameters_run | awk '{print $2}'`
set working_directory = `pwd`
set SCRATCH = ../scratch

#for full dataset, double last

@ full_last = 2 * $last

mainloop:

\rm ${model_a}_$start.par >& /dev/null
\rm ${model_b}_$start.par >& /dev/null

@ prev = $start - 1

if ($mode == 0) goto reconstruct

set firstn = $first
@ lastn = $first + $incr - 1

while ( $lastn <= $last )

  ${working_directory}/mrefine_n_halves_v9p07_helical.com $firstn $lastn $start $model_a $stack_a &

  if ( $lastn == $last ) then
    goto alignment_done_a
  endif
  @ firstn = $firstn + $incr
  @ lastn = $lastn + $incr
  if ( $firstn >= $last ) set firstn = $last
  if ( $lastn >= $last ) set lastn = $last
end

alignment_done_a:

  set firstn = $first
  @ lastn = $first + $incr - 1

checkdone_a:

    sleep 5
    while ( $firstn <= $last )
	
      grep --binary-files=text "overall score" $SCRATCH/${model_a}_${start}.par_${firstn}_$lastn >& /dev/null
      if ($status) goto checkdone_a

      echo "Half A particles $firstn to $lastn, finished....  "`date`
      if ($firstn == $first ) head -65 $SCRATCH/${model_a}_${start}.par_${firstn}_${lastn} | grep --binary-files=text C >> ${working_directory}/${model_a}_${start}.par

      grep -v C --binary-files=text $SCRATCH/${model_a}_${start}.par_${firstn}_${lastn} >> ${working_directory}/${model_a}_${start}.par
      \rm $SCRATCH/${model_a}_${start}.par_${firstn}_${lastn} >& /dev/null

      @ firstn = $firstn + $incr
      @ lastn = $lastn + $incr
      if ( $lastn >= $last ) set lastn = $last
    end

set firstn = $first
@ lastn = $first + $incr - 1

while ( $lastn <= $last )

  ${working_directory}/mrefine_n_halves_v9p07_helical.com $firstn $lastn $start $model_b $stack_b &

  if ( $lastn == $last ) then
    goto alignment_done_b
  endif
  @ firstn = $firstn + $incr
  @ lastn = $lastn + $incr
  if ( $firstn >= $last ) set firstn = $last
  if ( $lastn >= $last ) set lastn = $last
end

alignment_done_b:

  set firstn = $first
  @ lastn = $first + $incr - 1

checkdone_b:

    sleep 5
    while ( $firstn <= $last )
	
      grep --binary-files=text "overall score" $SCRATCH/${model_b}_${start}.par_${firstn}_$lastn >& /dev/null
      if ($status) goto checkdone_b

      echo "Half B particles $firstn to $lastn, finished....  "`date`
      if ($firstn == $first ) head -65 $SCRATCH/${model_b}_${start}.par_${firstn}_${lastn} | grep --binary-files=text C >> ${working_directory}/${model_b}_${start}.par

      grep -v C --binary-files=text $SCRATCH/${model_b}_${start}.par_${firstn}_${lastn} >> ${working_directory}/${model_b}_${start}.par
      \rm $SCRATCH/${model_b}_${start}.par_${firstn}_${lastn} >& /dev/null

      @ firstn = $firstn + $incr
      @ lastn = $lastn + $incr
      if ( $lastn >= $last ) set lastn = $last
    end

reconstruct:
echo "Calculating 3D structure Half A...."

  \rm $SCRATCH/${model_a}_${prev}.mrc*
  \rm $SCRATCH/${model_a}_${prev}.par*
  
  \rm $SCRATCH/${model_a}_mreconstruct.log
  \rm $SCRATCH/${model_a}.shft_* 
  \rm $SCRATCH/${model_a}_mrefine_n.log_*
  ${working_directory}/mreconstruct_halves_v9p07_helical.com $first $last $start $model_a $stack_a&

checkdoner_a:

  sleep 2

  grep --binary-files=text "mreconstruct.com finished" $SCRATCH/${model_a}_mreconstruct.log >& /dev/null
  if ($status) goto checkdoner_a

  ls ${working_directory}/${model_a}_${start}.mrc >& /dev/null
  if ($status) goto checkdoner_a

  sleep 1

  cat $SCRATCH/${model_a}_${start}.res >> ${working_directory}/${model_a}_${start}.par
  \rm $SCRATCH/${model_a}_${start}.res
  \rm $SCRATCH/${model_a}_${start}.par


echo "Calculating 3D structure Half B...."

  \rm $SCRATCH/${model_b}_${prev}.mrc*
  \rm $SCRATCH/${model_b}_${prev}.par*
  
  \rm $SCRATCH/${model_b}_mreconstruct.log
  \rm $SCRATCH/${model_b}.shft_* 
  \rm $SCRATCH/${model_b}_mrefine_n.log_*
  ${working_directory}/mreconstruct_halves_v9p07_helical.com $first $last $start $model_b $stack_b&

checkdoner_b:

  sleep 2

  grep --binary-files=text "mreconstruct.com finished" $SCRATCH/${model_b}_mreconstruct.log >& /dev/null
  if ($status) goto checkdoner_b

  ls ${working_directory}/${model_b}_${start}.mrc >& /dev/null
  if ($status) goto checkdoner_b

  sleep 1

  cat $SCRATCH/${model_b}_${start}.res >> ${working_directory}/${model_b}_${start}.par
  \rm $SCRATCH/${model_b}_${start}.res
  \rm $SCRATCH/${model_b}_${start}.par

#use CombineHalfParamsFre.py to get mixed parameter file for full dataset 
if ($mode == 0) then
		grep -v C ${model_a}_${prev}.par > half_a_stripped.par
		grep -v C ${model_b}_${prev}.par > half_b_stripped.par
		/home/galushin/goldhelix_sandbox/finished_scripts/CombineHalfParamsFre.py -a half_a_stripped.par -b half_b_stripped.par
		mv temp_halves_combined.par ${model_m}_${prev}.par
	else
		grep -v C ${model_a}_${start}.par > half_a_stripped.par
		grep -v C ${model_b}_${start}.par > half_b_stripped.par
		/home/galushin/goldhelix_sandbox/finished_scripts/CombineHalfParamsFre.py -a half_a_stripped.par -b half_b_stripped.par
		mv temp_halves_combined.par ${model_m}_${start}.par
	endif

\rm half_a_stripped.par
\rm half_b_stripped.par


echo "Calculating 3D structure Combined Halves...."

  \rm $SCRATCH/${model_m}_${prev}.mrc*
  \rm $SCRATCH/${model_m}_${prev}.par*
  
  \rm $SCRATCH/${model_m}_mreconstruct.log
  \rm $SCRATCH/${model_m}.shft_* 
  \rm $SCRATCH/${model_m}_mrefine_n.log_*
  ${working_directory}/mreconstruct_halves_v9p07_helical.com $first $full_last $start $model_m $stack_m&

checkdoner_m:

  sleep 2

  grep --binary-files=text "mreconstruct.com finished" $SCRATCH/${model_m}_mreconstruct.log >& /dev/null
  if ($status) goto checkdoner_m

  ls ${working_directory}/${model_m}_${start}.mrc >& /dev/null
  if ($status) goto checkdoner_m

  sleep 1

  cat $SCRATCH/${model_m}_${start}.res >> ${working_directory}/${model_m}_${start}.par
  \rm $SCRATCH/${model_m}_${start}.res
  \rm $SCRATCH/${model_m}_${start}.par

sleep 10

echo "Cycle ${start} finished.... "`date`

if ($start < $end ) then
  @ start = $start + 1
  goto mainloop
endif
date
