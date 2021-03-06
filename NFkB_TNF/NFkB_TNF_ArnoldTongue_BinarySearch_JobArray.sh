#!/bin/bash

source init_env.sh

index=$1
output=$2
volume=$3
input_periods=$4
tolerance=$5
Ntrials=$6
population_average=$7

echo "index=$index, output=$output, volume=$volume, input_periods=$input_periods, tolerance=$tolerance, Ntrials=$Ntrials, population_average=$population_average"

matlab -nodisplay -singleCompThread -r "NFkB_TNF_ArnoldTongue_BinarySearch_JobArrayWrapper($index, '$output', $volume, $input_periods, $tolerance, $Ntrials, $population_average); exit;"
