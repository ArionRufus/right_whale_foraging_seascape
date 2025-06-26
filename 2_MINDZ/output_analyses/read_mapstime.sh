#!/bin/bash
# Argument passed to the script
#SBATCH --job-name=time_maps                 # nom de la tache à runer
#SBATCH --output=time_maps.out               # nom du fichier dans lequel on peut lire les infos d'outputs
#SBATCH --partition=batch                    # nom de la partition a choisir (sinfo)
#SBATCH --nodes=1                            # node count
#SBATCH --cpus-per-task=1                    # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=64G                    # total memory per node (4 GB per cpu-core is default)
#SBATCH --time=00-07:30:00                   # total run time limit (DD-HH:MM:SS)
#SBATCH --account=def-frmap5
#SBATCH --mail-type=ALL                      # send an email in all cases (job started, job ended, job aborted)
#SBATCH --mail-user=andeol.bourgouin.1@ulaval.ca

# load required modules:
module load  StdEnv/2020  gcc/9.3.0  udunits/2.2.28  gdal/3.5.1  r/4.2.2

#installer packages
Rscript -e 'renv::load()'  #pour aller chercher renv.lock
Rscript -e "renv::activate()"

# executer le script:
Rscript time_maps.R
