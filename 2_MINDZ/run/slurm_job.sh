#!/bin/bash 
#SBATCH --account=def-fmaps
#SBATCH -t 00-8:00  
#SBATCH --mem-per-cpu=16G  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cebrennan.climate@gmail.com

PROGNOM=`basename $0`

usage()
{
  echo usage:  "sbatch {PROGNOM} executable.x"
  echo ' '
}

# Check for # of tokens

if [ $# -lt 1 ]
then
  usage
  exit 1
elif [$# -eq 2]
then
  echo "${PROGRNOM}=run the MINDZ (Model of INDividualZooplankton} executable."
  echo ''
  export OMP_NUM_THREADS=$2
  ./$1
else 
  echo "${PROGRNOM}=run the MINDZ (Model of INDividualZooplankton} executable"
  echo ''
  ./$1
fi
