#!/bin/bash

module load matlab

index=$1
output_prefix=$2

echo "index=$index, output_prefix=$output_prefix"

matlab -nodisplay -singleCompThread -r "VanDerPol_ArnoldTongue_Stochastic_JobArray($index, $output_prefix)"
