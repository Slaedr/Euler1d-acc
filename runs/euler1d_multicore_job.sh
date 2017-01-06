#! /bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=00:05:00
#PBS -A rck-371-aa
#PBS -o sodst_acc_cpus-2threads.log
#PBS -e sodst_acc_cpus.err
#PBS -N sodst_acc_cpus

# to be executed from Euler1d-acc/runs/ or equivalent

module load CUDA/7.5.18
module load PGI/16.3-GCC-4.9.3-2.25
#export PGI_ACC_TIME=1
export ACC_NUM_CORES=2

cd $PBS_O_WORKDIR
../build/euler1d-acc ../cases/sodst.control
