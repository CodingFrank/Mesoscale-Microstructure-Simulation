#!/bin/bash
#
# srun command for CSCI-6360 group project.
# USAGE: /full/path/to/./q_MC.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --mail-type=END
#SBATCH --mail-user=zhengc2@rpi.edu
#SBATCH -D /gpfs/u/scratch/ACME/ACMEchzh/PPC/PROJ
#SBATCH --partition debug
#SBATCH -t 20
#SBATCH -N 8
#SBATCH -n 256
#SBATCH --overcommit
#SBATCH -o /gpfs/u/scratch/ACME/ACMEchzh/PPC/PROJ/output1.log


# Use mc.out <Number of threads per rank> <Domain size> <Number of updates> <Output Interval> <Initial grian size>
srun --runjob-opts="--mapping TEDCBA" /gpfs/u/scratch/ACME/ACMEchzh/PPC/PROJ/mc.out 0 1024 1 10 20
