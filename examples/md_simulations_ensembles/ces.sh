. ../set_environment.sh
ln -s $ENCORE_BUILD .
ln -s $MDANALYSIS_BUILD .

# Trajectories containing the ensembles to be compared
trajs=(traj_c22-star.xtc  traj_c36.xtc  traj_ff99sb-ildn-star.xtc)

# number of ensembles
nensembles=3

# command line. Options:
# --topology: topology for the provided trajectories
# --mode=ces: clustering ensemble similarity
# --nensembles: numbers of ensembles that will be used for calculation
# --save-matrix: save the similarity (-RMSD) matrix to file
# --superimpose: best superimpose structures before calculating the RMSD
# --convergence, --maxiter, --preferences, --lambda: parameters for the Affinity propagation clustering algorithm (see help)
# --np: maximum number of cores to be used 


cmdline="./similarity.py --save-matrix=minusrmsd_pw.npz --superimpose --convergence=50 --maxiter=500 --preferences=-5.0 --np=$NP --nensembles $nensembles --topology topology.pdb --mode=ces -v --lambda=0.9"

# add the --ensemble*-trajectory options (--ensemble1-trajectory traj_c22-star.xtc --ensemble2-trajectory traj_c36.xtc etc.)
for i in $(seq 0 2); do
	cmdline=$cmdline" --ensemble$((i+1))-trajectory ${trajs[$i]}"
done

echo "Now running: $cmdline"
echo  $ENCORE_BUILD/similarity.py
$PYTHON encore/similarity.py $cmdline 

