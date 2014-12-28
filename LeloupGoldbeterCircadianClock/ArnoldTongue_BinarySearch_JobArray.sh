#!/bin/bash

source init_env.sh

index=$1
output=$2
input_periods=$3
tolerance=$4
Ntrials=$5
population_average=$6
entrainment_threshold=$7
max_input_amplitude=$8
final_time=$9
transient_time_factor=$10
num_of_timepoints=$11

echo "index=$index, output=$output, input_periods=$input_periods, tolerance=$tolerance, Ntrials=$Ntrials, population_average=$population_average, entrainment_threshold=$entrainment_threshold, max_input_amplitude=$max_input_amplitude, final_time=$final_time, transient_time_factor=$transient_time_factor, num_of_timepoints=$num_of_timepoints"

matlab -nodisplay -singleCompThread -r "ArnoldTongue_BinarySearch_JobArray($index, '$output', $input_periods, $tolerance, $Ntrials, $population_average, $entrainment_threshold, $max_input_amplitude, $final_time, $transient_time_factor, $num_of_timepoints); exit;"
