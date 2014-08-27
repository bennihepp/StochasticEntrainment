#!/bin/bash

module load matlab

index=$1
output=$2

echo "index=$index, output=$output"

matlab -nodisplay -singleCompThread -r "VanDerPol_ArnoldTongue_Stochastic_BinarySearch_JobArray($index, '$output'); exit;"
