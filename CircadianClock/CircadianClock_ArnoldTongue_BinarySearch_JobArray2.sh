#!/bin/bash

source init_env.sh

index=$1
output=$2
input_periods=$3
tolerance=$4
Ntrials=$5
population_average=$6
Nthreads=$7

echo "index=$index, output=$output, input_periods=$input_periods, tolerance=$tolerance, Ntrials=$Ntrials, population_average=$population_average, Nthreads=$Nthreads"

matlab -nodisplay -singleCompThread -r "CircadianClock_ArnoldTongue_BinarySearch_JobArray2Wrapper($index, '$output', $input_periods, $tolerance, $Ntrials, $population_average, $Nthreads); exit;"
