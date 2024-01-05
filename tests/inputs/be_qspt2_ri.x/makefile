clean:
	rm -f log pbs.* RESTART ENERGY_QP *.dat *.out

profile:
	nvprof                   /home/ymbyun/programs/molgw-1.F-2019.03.20-N17.O1-GW.inner/molgw molgw.in > log
#	nvprof --print-gpu-trace /home/ymbyun/programs/molgw-1.F-2019.03.20-N17.O1-GW.inner/molgw molgw.in > log
#	nvprof --print-api-trace /home/ymbyun/programs/molgw-1.F-2019.03.20-N17.O1-GW.inner/molgw molgw.in > log

profile_visual:
	nvprof --analysis-metrics -o be_qspt2_analysis.nvvp /home/ymbyun/programs/molgw-1.F-2019.03.20-N17.O1-GW.inner/molgw molgw.in > log
#	nvprof                    -o be_qspt2.nvvp          /home/ymbyun/programs/molgw-1.F-2019.03.20-N17.O1-GW.inner/molgw molgw.in > log

run:
	/home/ymbyun/programs/molgw-1.F-2019.03.20-N17.O1-GW.inner/molgw molgw.in > log 2>&1
