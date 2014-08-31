#!/bin/bash

source init_env.sh

index=$1
output=$2
input_periods=$3
tolerance=$4
population_average=$5

echo "index=$index, output=$output, input_periods=$input_periods, tolerance=$tolerance, population_average=$population_average"

matlab -nodisplay -singleCompThread -r "CircadianClock_ArnoldTongue_BinarySearch_JobArrayWrapper($index, '$output', $input_periods, $tolerance, $population_average); exit;"
