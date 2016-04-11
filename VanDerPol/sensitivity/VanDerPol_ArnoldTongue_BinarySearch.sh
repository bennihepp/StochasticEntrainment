#!/bin/bash

source init_env.sh

output=$1
volume=$2
input_periods=$3
tolerance=$4
Ntrials=$5
population_average=$6
entrainment_threshold=$7
parameterB=$8
parameterD=$9

echo "output=$output, volume=$volume, input_periods=$input_periods, tolerance=$tolerance, Ntrials=$Ntrials, population_average=$population_average, entrainment_threshold=$entrainment_threshold, parameterB=$parameterB, parameterD=$parameterD"

matlab -nodisplay -singleCompThread -r "VanDerPol_ArnoldTongue_BinarySearch('$output', $volume, $input_periods, $tolerance, $Ntrials, $population_average, $entrainment_threshold, $parameterB, $parameterD); exit;"
