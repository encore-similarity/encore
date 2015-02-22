. ../set_environment.sh
ln -s ../$ENCORE_BUILD .
ln -s ../$MDANALYSIS_BUILD .

# List of files containing the ensembles to be compared. The x-ray derived ensemble is included as a .DCD trajectory file.
pdbs=(filtered_1XQQ.pdb  filtered_2K39.pdb  filtered_2KOX.pdb filtered_2LJ5.pdb filtered_2NR2.pdb xray.dcd)

# command line
nensembles=6

# Command line. Please notice the --save-matrix line, which saves on disk the -RMSD matrix that we use as a similarity matrix for affinity propagation. We will re-use it later in the dimred.sh.
cmdline="./similarity.py --save-matrix=minusrmsd_pw.npz --superimpose --convergence=50 --maxiter=500 --preferences=-4.0 --np=$NP --nensembles $nensembles --topology filtered_reference.pdb --mode=clustering -v"

# for each ensemble: add the --ensemble-trajectory option to the command line. For instance: --ensemble1-trajectory filtered_1XQQ.pdb --ensemble2-trajectory filtered_2K39.pdb etc.
for i in $(seq 0 5); do
	cmdline=$cmdline" --ensemble$((i+1))-trajectory ${pdbs[$i]}"
done

# run encore
echo "Now running: $cmdline"
echo "Results will be written in file: clustering.log"
echo  $ENCORE_BUILD/similarity.py
$PYTHON encore/similarity.py $cmdline &> clustering.log

