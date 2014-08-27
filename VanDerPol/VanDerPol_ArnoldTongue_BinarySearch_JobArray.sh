#!/bin/bash

module load matlab

index=$1
output=$2
population_average=$3

echo "index=$index, output=$output, population_average=$population_average"

matlab -nodisplay -singleCompThread -r "VanDerPol_ArnoldTongue_BinarySearch_JobArray($index, '$output', $population_average); exit;"
