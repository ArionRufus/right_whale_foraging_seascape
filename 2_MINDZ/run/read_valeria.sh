#!/bin/bash
# Argument passed to the script
#SBATCH --job-name=MINDZ                     # nom de la tache à runer
#SBATCH --output=MINDZ_execution.out         # nom du fichier dans lequel on peut lire les infos d'outputs
#SBATCH --partition=batch_96h                    # nom de la partition a choisir (sinfo)
#SBATCH --nodes=1                            # node count
#SBATCH --cpus-per-task=1                    # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=64G                    # total memory per node (4 GB per cpu-core is default)
#SBATCH --time=03-00:00:00                   # total run time limit (DD-HH:MM:SS)
#SBATCH --account=def-frmap5
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=ALL                      # send an email in all cases (job started, job ended, job aborted)
#SBATCH --mail-user=andeol.bourgouin.1@ulaval.ca

# load required modules:
module load netcdf-fortran python

cd MINDZ_andeol/MINDZ/run

# executer le script:
./multiple_dates.py

