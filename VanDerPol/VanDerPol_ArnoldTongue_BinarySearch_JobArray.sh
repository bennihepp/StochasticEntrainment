#!/bin/bash

source init_env.sh

index=$1
output=$2
input_periods=$3
tolerance=$4
Ntrials=$5
population_average=$6
entrainment_threshold=$7

echo "index=$index, output=$output, input_periods=$input_periods, tolerance=$tolerance, Ntrials=$Ntrials, population_average=$population_average, entrainment_threshold=$entrainment_threshold"

matlab -nodisplay -singleCompThread -r "VanDerPol_ArnoldTongue_BinarySearch_JobArray($index, '$output', $input_periods, $tolerance, $Ntrials, $population_average, $entrainment_threshold); exit;"
