. ../set_environment.sh
ln -s $ENCORE_BUILD .
ln -s $MDANALYSIS_BUILD .

# Just one ensemble, since we want to evaluate convergence.                                                                                         
nensembles=1

# command line. Options:
# --topology: topology for the provided trajectories
# --mode=dres: dimensionality reduction ensemble similarity
# --nensembles: numbers of ensembles that will be used for calculation
# --load-matrix: load similarity (-RMSD) matrix instead of calculating it
# --change-matrix-sign: change the sign of the elements of the loaded matrix (so that we have a +RMSD matrix)
# --use-density, --dim, --samples, --nsteps=, --ncycle, --neighborhood-cutoff, --spe-mode: parameters for the Stochastic proximity embedding algorithm (see help for details)
# 
# --evaluate-convergence: do evaluate convergence on ensemble 1
# --window-size: number of frames to be used for the first window. The second and so on window will be multiples of the first one. (i.e. with --window-size=25 windows will be long 25, 50, 75, 100 ... frames)

cmdline="./similarity.py --load-matrix=minusrmsd_ff99sb-ildn-star.npz --change-matrix-sign --superimpose --mode=dres --dim=2 --samples=10000 --nsteps=1000000 --ncycle=100 --neighborhood-cutoff=1.7 --spe-mode=vanilla --np=$NP --nensembles $nensembles --topology topology.pdb --evaluate-convergence --evaluate-convergence-mode=increasing-window --window-size=25 -v"

# Run the same command line for both trajectories
for ff in ff99sb-ildn-star c22-star; do
	traj=traj_"$ff".xtc

	cmdline=$cmdline" --ensemble1-trajectory="$traj

	echo "Now running: $ff"
	echo  $ENCORE_BUILD/similarity.py
	$PYTHON encore/similarity.py $cmdline
done
