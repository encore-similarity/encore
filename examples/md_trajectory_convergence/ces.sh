. ../set_environment.sh
ln -s $ENCORE_BUILD .
ln -s $MDANALYSIS_BUILD .

# Just one ensemble, since we want to evaluate convergence.
nensembles=1


# Command line. Options:
# --topology: topology for the provided trajectories
# --mode=ces: clustering ensemble similarity
# --nensembles: numbers of ensembles that will be used for calculation
# --save-matrix: save the similarity (-RMSD) matrix to file
# --superimpose: best superimpose structures before calculating the RMSD
# --convergence, --maxiter, --preferences, --lambda: parameters for the Affinity              propagation clustering algorithm (see help)
# --np: maximum number of cores to be used
#
# --evaluate-convergence: do evaluate convergence on ensemble 1
# --window-size: number of frames to be used for the first window. The second and so on window will be multiples of the first one. (i.e. with --window-size=25 windows will be long 25, 50, 75, 100 ... frames)

cmdline="./similarity.py --save-matrix=minusrmsd_ff99sb-ildn-star.npz --superimpose --convergence=50 --maxiter=500 --preferences=-2 --np=$NP --nensembles $nensembles --topology topology.pdb --mode=ces -v --evaluate-convergence --evaluate-convergence-mode=increasing-window --window-size=25"

# Run the same command line for both trajectories
for ff in ff99sb-ildn-star c22-star; do
	traj=traj_"$ff".xtc

	cmdline=$cmdline" --ensemble1-trajectory="$traj

	echo "Now running: $ff"
	echo  $ENCORE_BUILD/similarity.py
	$PYTHON encore/similarity.py $cmdline
done
