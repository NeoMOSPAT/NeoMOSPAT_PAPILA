#!/bin/bash

#####################################
folder="/home/m/m301017/Codes-PAPILA/NeoMOSPAT/SIMULDIR/TS_models_20230704/"
folder=/home/m/m301017/Codes-PAPILA/NeoMOSPAT/SIMULDIR/New_Images/
folder="/home/m/m301017/Codes-PAPILA/NeoMOSPAT/SIMULDIR/2024/"
folder="/home/m/m301017/Codes-PAPILA/NeoMOSPAT/SIMULDIR/TS_ensembles_20230704/"
#TS_ensembles-MAPLA_20230704/
dependency="--dependency=after:8825057+1"
#dependency=""

for file in $folder/*.py
do
    echo "Submitting file $file"
    ## Submit job and retrieve jobidi
    jobname=$(basename $file)
    sentence=$(sbatch $dependency -J $jobname job.sh $file)
    #echo sbatch $dependency -J $jobname job.sh $file
    stringarray=($sentence)
    jobid=(${stringarray[3]})

    #echo "Submitting file $file"
    echo "sbatch $dependency -J $jobname job.sh $file"
    #echo $jobid
    #echo $sentence
    #echo " "
    ## Create dependecy line based on jobid for next file in recurrsion
    dependency="--dependency=after:$jobid+1"
    
done
