#!/bin/bash

source init_env.sh

index=$1
output=$2
input_periods=$3
population_average=$4

echo "index=$index, output=$output, input_periods=$input_periods, population_average=$population_average"

matlab -nodisplay -singleCompThread -r "CircadianClock_ArnoldTongue_BinarySearch_JobArrayWrapper($index, '$output', $input_periods, $population_average); exit;"
