#!/bin/bash
#$ -V                    # export all environment variables
#$ -S /bin/bash          # command interpreter to be used
#$ -N layer1234           # job name
#$ -q all.q              # queue type
#$ -pe mpi_48 48         # request slot range for parallel jobs mpi_(core per node) (total core)
#$ -j Y                  # merge e and o file and pe and po file
#$ -o $JOB_NAME.o$JOB_ID # merge po and o file
#$ -cwd                  # use current working directory
###$ -t 1-10:4 # 1, 5, 9 # job array. It can be reached using $SGE_TASK_ID 

###### don't touch below #####
echo "Got $NSLOTS slots."
cat $TMPDIR/machines
export OMP_NUM_THREADS=1
##############################

# go to working directory
cd $SGE_O_WORKDIR
VASP="/home/shared/programs/vasp/5.4.4+vtst+vaspsol+beef/vasp.5.4.4/bin/vasp_std"
#VASP="/home/shared/programs/vasp/5.4.4+vtst+vaspsol+beef/vasp.5.4.4/bin/vasp_gam"
mpirun -machinefile $TMPDIR/machines -n $NSLOTS $VASP
# keep machinefile option above
