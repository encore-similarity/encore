. ../set_environment.sh
ln -s $ENCORE_BUILD .
ln -s $MDANALYSIS_BUILD .

# Just one ensemble, since we want to evaluate convergence.                                                                                         
nensembles=1

# Command line
cmdline="./similarity.py --load-matrix=minusrmsd_ff99sb-ildn-star.npz --change-matrix-sign --superimpose --mode=dimred --dim=2 --samples=10000 --nsteps=1000000 --ncycle=100 --neighborhood-cutoff=1.7 --spe-mode=vanilla --np=$NP --nensembles $nensembles --topology topology.pdb --evaluate-convergence --evaluate-convergence-mode=increasing-window --window-size=100 -v"

# Run the same command line for both trajectories
for ff in ff99sb-ildn-star c22-star; do
	traj=traj_"$ff".xtc

	cmdline=$cmdline" --ensemble1-trajectory="$traj

	echo "Now running: $ff"
	echo "Results will be written in file: dimred_"$ff".log"
	echo  $ENCORE_BUILD/similarity.py
	$PYTHON encore/similarity.py $cmdline
done
