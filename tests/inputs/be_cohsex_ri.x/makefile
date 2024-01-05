clean:
	rm -f log pbs.* RESTART ENERGY_QP *.dat *.out

run:
	/home/ymbyun/programs/molgw-1.F-2019.03.20-N17.O1/molgw molgw.in > log 2>&1

#export MKL_NUM_THREADS=4
#export OMP_NUM_THREADS=4
#/home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.GW.Z-4.AE.5Z-5.S-6.OpenMP-N17.O1/molgw molgw.in >& log

#mpiexec -n 1 /home/ybyun/programs/molgw-1.F-1.HF.alpha-2.Vxc-3.GW.Z-4.AE.5Z-5.S-N17.O1/molgw molgw.in >& log
