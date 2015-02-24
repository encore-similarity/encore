. ../set_environment.sh
ln -s $ENCORE_BUILD .
ln -s $MDANALYSIS_BUILD .

# Trajectories containing the ensembles to be compared
trajs=(traj_c22-star.xtc  traj_c36.xtc  traj_ff99sb-ildn-star.xtc)

# number of ensembles
nensembles=3

# command line
cmdline="./similarity.py --save-matrix=minusrmsd_pw.npz --superimpose --convergence=50 --maxiter=500 --preferences=-5.0 --np=$NP --nensembles $nensembles --topology topology.pdb --mode=clustering -v --lambda=0.9"

# add the --ensemble*-trajectory options (--ensemble1-trajectory traj_c22-star.xtc --ensemble2-trajectory traj_c36.xtc etc.)
for i in $(seq 0 2); do
	cmdline=$cmdline" --ensemble$((i+1))-trajectory ${trajs[$i]}"
done

echo "Now running: $cmdline"
echo "Results will be written in file: clustering.log"
echo  $ENCORE_BUILD/similarity.py
$PYTHON encore/similarity.py $cmdline 

