. ../set_environment.sh
ln -s $ENCORE_BUILD .
ln -s $MDANALYSIS_BUILD .

# List of files containing the ensembles to be compared. The x-ray derived ensemble is included as a .DCD trajectory file.
pdbs=(filtered_1XQQ.pdb  filtered_2K39.pdb  filtered_2KOX.pdb filtered_2LJ5.pdb filtered_2NR2.pdb xray.dcd)

# number of ensembles
nensembles=6


# command line. Options:
# --topology: topology for the provided trajectories
# --mode=ces: clustering ensemble similarity
# --nensembles: numbers of ensembles that will be used for calculation
# --save-matrix: save the similarity (-RMSD) matrix to file
# --superimpose: best superimpose structures before calculating the RMSD
# --convergence, --maxiter, --preferences, --lambda: parameters for the Affinity              propagation clustering algorithm (see help)
# --np: maximum number of cores to be used

cmdline="./similarity.py --save-matrix=minusrmsd_pw.npz --superimpose --convergence=50 --maxiter=500 --preferences=-5.0 --np=$NP --nensembles $nensembles --topology filtered_reference.pdb --mode=ces -v"

# for each ensemble: add the --ensemble-trajectory option to the command line. For instance: --ensemble1-trajectory filtered_1XQQ.pdb --ensemble2-trajectory filtered_2K39.pdb etc.
for i in $(seq 0 5); do
	cmdline=$cmdline" --ensemble$((i+1))-trajectory ${pdbs[$i]}"
done

# run encore
echo "Now running: $cmdline"
echo  $ENCORE_BUILD/similarity.py
$PYTHON encore/similarity.py $cmdline

