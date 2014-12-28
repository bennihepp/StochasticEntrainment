#!/bin/bash

OUTPUT_FOLDER=output_50/
JOBID=lg
POPULATION_AVERAGE=false
#INPUT_PERIODS=20:0.1:28
#NUM_OF_PERIODS=81
INPUT_PERIODS=18:0.5:30
NUM_OF_PERIODS=25
AMPLITUDE_TOLERANCE=1e-2
NTRIALS=50
ENTRAINMENT_THRESHOLD=0.9
MAX_INPUT_AMPLITUDE=0.3
FINAL_TIME=$((200*24))
TRANSIENT_TIME_FACTOR=0.1
NUM_OF_TIMESTEPS=50001
export NUM_OF_THREADS=16
NUM_OF_CPUS=${NUM_OF_THREADS}
RUNTIME="12:00"
#RUNTIME="1:00"

mkdir $OUTPUT_FOLDER
mkdir $OUTPUT_FOLDER/logs

#bsub -n 1 -R "rusage[mem=1024]" -W 1:00 -J "${JOBID}a" -o $OUTPUT_FOLDER/logs/VanDerPol_ArnoldTongue_BinarySearch_JobArray2_Init.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray2.sh -1 $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE
bash ArnoldTongue_BinarySearch_JobArray.sh -1 $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE $ENTRAINMENT_THRESHOLD $MAX_INPUT_AMPLITUDE $FINAL_TIME $TRANSIENT_TIME_FACTOR $NUM_OF_TIMESTEPS

#bsub -n ${NUM_OF_CPUS} -R "rusage[mem=2048]" -W ${RUNTIME} -J "${JOBID}b[1-$NUM_OF_PERIODS]" -o $OUTPUT_FOLDER/logs/ArnoldTongue_BinarySearch_JobArray_%I.out bash ArnoldTongue_BinarySearch_JobArray.sh "\$LSB_JOBINDEX" $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE $ENTRAINMENT_THRESHOLD $MAX_INPUT_AMPLITUDE $FINAL_TIME $TRANSIENT_TIME_FACTOR $NUM_OF_TIMESTEPS
#JOBINDEX=38; bsub -n ${NUM_OF_CPUS} -R "rusage[mem=2048]" -W ${RUNTIME} -J "${JOBID}b${JOBINDEX}" -o $OUTPUT_FOLDER/logs/VanDerPol_ArnoldTongue_BinarySearch_JobArray2_${JOBINDEX}.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray2.sh ${JOBINDEX} $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE

#bsub -n 1 -R "rusage[mem=1536]" -W 1:00 -w "done(${JOBID}b)" -J "${JOBID}c" -o $OUTPUT_FOLDER/logs/VanDerPol_ArnoldTongue_BinarySearch_JobArray2_Combine.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray2.sh 0 $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE
#bash ArnoldTongue_BinarySearch_JobArray.sh 0 $OUTPUT_FOLDER "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE $ENTRAINMENT_THRESHOLD $MAX_INPUT_AMPLITUDE $FINAL_TIME $TRANSIENT_TIME_FACTOR $NUM_OF_TIMESTEPS

