#!/bin/bash
#SBATCH -J Neomospat
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --mem=0
#SBATCH --output=logs/log_%j.out
#SBATCH --account=bm1234
#SBATCH --mail-type=FAIL
#SBATCH --time=05:00:00
set -e
ulimit -s 204800
ulimit -c 0

ml python3

## If job is passed without an argument then it searches for IncludeFile.py 
python -W ignore -u main.py $1 #2>&1 | tee log.txt

