#!/bin/bash -l

#PBS -q login
#PBS -l nodes=1:ppn=8
#PBS -l mem=48GB
#PBS -l walltime=10:00:00
#PBS -N MP2.ri.o
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
#/home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-N17.O1-bug.fix/molgw_omp_aff.x molgw.in > log_omp
#/home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-8.COHSEX.QSGW-N17.O1/molgw_omp molgw.in > log
#/home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-8.COHSEX.QSGW-N17.O1-2019.removed/molgw_omp molgw.in > log
#/home/ybyun/programs/molgw-2.A-2018.11.24-N17.O1/molgw_omp_aff.x_0 molgw.in > log
#/home/ybyun/programs/molgw-2.A-2019.03.19-N17.O1/molgw_omp molgw.in > log
#mpiexec /home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-N17.O1-bug.fix/molgw_mpi_ymb.o molgw.in > log_mpi_ymb.ooo
#mpiexec /home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-N17.O1-bug.fix/molgw_mpi_ymb.x molgw.in > log_mpi_ymb.xxx
#mpiexec /home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.ERI.i8-4.S-5.OpenMP-6.RPA-7.Eqp.A.Z-8.COHSEX.QSGW-N17.O1/molgw_mpi molgw.in > log
mpiexec /home/ybyun/programs/molgw-1.F-2019.03.20-N17.O1/molgw_mpi molgw.in > log
date

echo ""
echo "Done!"
