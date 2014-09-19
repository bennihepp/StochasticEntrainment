#!/bin/bash

source init_env.sh

index=$1
output=$2
volume=$3
input_periods=$4
tolerance=$5
Ntrials=$6
population_average=$7
entrainment_threshold=$8

echo "index=$index, output=$output, volume=$volume, input_periods=$input_periods, tolerance=$tolerance, Ntrials=$Ntrials, population_average=$population_average, entrainment_threshold=$entrainment_threshold"

matlab -nodisplay -singleCompThread -r "NFkB_TNF_ArnoldTongue_BinarySearch_JobArray2Wrapper($index, '$output', $volume, $input_periods, $tolerance, $Ntrials, $population_average, $entrainment_threshold); exit;"
