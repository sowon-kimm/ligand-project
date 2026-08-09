#!/bin/sh
#PBS -N mag_add
#PBS -V
#PBS -q normal
#PBS -A vasp
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -l walltime=48:00:00
#PBS -m abe 
#PBS -M ksone23@kentech.ac.kr

###### don't touch below #####
#echo "Got $NSLOTS slots."
#cat $TMPDIR/machines
#export OMP_NUM_THREADS=1
##############################

# go to working directory
#cd $SGE_O_WORKDIR
#VASP="/home01/x3370a02/vasp_std"

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

module purge
module load craype-mic-knl
module load intel impi

VASP="/home01/x3370a02/vasp_std"
mpirun $VASP > vasp.out

#module purge
#module load craype-mic-knl
#./test.exe
