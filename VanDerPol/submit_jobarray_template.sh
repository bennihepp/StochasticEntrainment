#!/bin/bash

# define some parameters for the job array
JOBID=vp1
export NUM_OF_THREADS=16
NUM_OF_CPUS=${NUM_OF_THREADS}
RUNTIME="12:00"

# define some input arguments for the VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh script
OUTPUT_FOLDER=output_500/
INPUT_PERIODS=4:0.2:18
NUM_OF_PERIODS=71
AMPLITUDE_TOLERANCE=1e-2
NTRIALS=500
POPULATION_AVERAGE=false
ENTRAINMENT_THRESHOLD=0.9

# create folder for results and log files
mkdir $OUTPUT_FOLDER
mkdir $OUTPUT_FOLDER/logs

# create info.mat file containing
bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh -1 $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE $ENTRAINMENT_THRESHOLD

# submit job array
#bsub -n ${NUM_OF_CPUS} -R "rusage[mem=2048]" -W ${RUNTIME} -J "${JOBID}[1-$NUM_OF_PERIODS]" -o $OUTPUT_FOLDER/logs/ArnoldTongue_BinarySearch_JobArray_%I.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh "\$LSB_JOBINDEX" $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE $ENTRAINMENT_THRESHOLD

# combine results into a single file: VanDerPol_ArnoldTongue_BinarySearch_JobArray_[...].mat
#bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh 0 $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE $ENTRAINMENT_THRESHOLD

