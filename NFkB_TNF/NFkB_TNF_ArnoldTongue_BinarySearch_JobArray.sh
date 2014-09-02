#!/bin/bash

source init_env.sh

index=$1
output=$2
input_periods=$3
tolerance=$4
Ntrials=$5
population_average=$6

echo "index=$index, output=$output, input_periods=$input_periods, tolerance=$tolerance, Ntrials=$Ntrials, population_average=$population_average"

matlab -nodisplay -singleCompThread -r "NFkB_TNF_ArnoldTongue_BinarySearch_JobArrayWrapper($index, '$output', $input_periods, $tolerance, $Ntrials, $population_average); exit;"
