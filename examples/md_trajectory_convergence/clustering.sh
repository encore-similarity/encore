. ../set_environment.sh
ln -s ../$ENCORE_BUILD .
ln -s ../$MDANALYSIS_BUILD .

# Just one ensemble, since we want to evaluate convergence.
nensembles=1

# Command line
cmdline="./similarity.py --save-matrix=minusrmsd_ff99sb-ildn-star.npz --superimpose --convergence=50 --maxiter=500 --preferences=-2 --np=$NP --nensembles $nensembles --topology topology.pdb --mode=clustering -v --evaluate-convergence --evaluate-convergence-mode=increasing-window --window-size=100"

# Run the same command line for both trajectories
for ff in ff99sb-ildn-star c22-star; do
	traj=traj_"$ff".xtc

	cmdline=$cmdline" --ensemble1-trajectory="$traj

	echo "Now running: $ff"
	echo "Results will be written in file: clustering_"$ff".log"
	echo  $ENCORE_BUILD/similarity.py
	$PYTHON encore/similarity.py $cmdline &> clustering_"$ff".log
done
