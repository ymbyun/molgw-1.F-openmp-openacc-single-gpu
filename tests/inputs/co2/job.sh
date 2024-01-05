#!/bin/bash -l

#PBS -q head
#PBS -l nodes=1:ppn=1
#PBS -l mem=24GB
#PBS -l walltime=500:00:00
#PBS -N bug.fix
#PBS -e pbs.err
#PBS -o pbs.out

cd $PBS_O_WORKDIR

# Long info
echo ""
env | grep PBS

# Short info
NPROCS=$(wc -l < $PBS_NODEFILE)
NNODES=$(uniq $PBS_NODEFILE | wc -l)
echo ""
printf "%s %2d \n" "Number of procs:" $NPROCS
printf "%s %2d \n" "Number of nodes:" $NNODES
echo ""
cat $PBS_NODEFILE
echo ""

# OpenMP
#export MKL_NUM_THREADS=8
#export OMP_NUM_THREADS=8
#export OMP_STACKSIZE=24M
#export OMP_DISPLAY_ENV=true
#export KMP_AFFINITY=verbose,none
#export KMP_AFFINITY=verbose,scatter
#export KMP_AFFINITY=verbose,compact
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread

# Run
date
mpiexec /home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-N17.O1-bug.fix/molgw_mpi_ymb.x molgw.in > log_ymb.xxx
#mpiexec /home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-N17.O1-bug.fix/molgw_mpi_ymb.o molgw.in > log_ymb.ooo
date

echo ""
echo "Done!"
