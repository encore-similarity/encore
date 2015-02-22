. ../set_environment.sh
ln -s ../$ENCORE_BUILD .
ln -s ../$MDANALYSIS_BUILD .

if [[ ! -e minusrmsd_pw.npz ]]; then
	echo "-RMSD Matrix file ( minusrmsd_pw.npz ) not found! Please run the clustering example first."
	exit
fi

# List of files containing the ensembles to be compared. The x-ray derived ensemble is included as a .DCD trajectory file.
pdbs=(filtered_1XQQ.pdb filtered_2K39.pdb filtered_2KOX.pdb filtered_2LJ5.pdb filtered_2NR2.pdb xray.dcd)

# number of ensembles
nensembles=6

# Command line to be used. Notice the --load-matrix option, which means that the -RMSD matrix is loaded from disk. --change-matrix-sign inverts the sign of the loaded matrix, since we need a distance matrix for Stochastic proximity embedding.
cmdline="./similarity.py --mode=dimred --superimpose --np=6 --load-matrix=minusrmsd_pw.npz --change-matrix-sign --nensembles $nensembles --topology filtered_reference.pdb --use-density=resample --dim=3 --samples=10000 --nsteps=1000000 --ncycle=100 --neighborhood-cutoff=1.5 --spe-mode=vanilla -v"

# for each ensemble: add the --ensemble-trajectory option to the command line. For instance: --ensemble1-trajectory filtered_1XQQ.pdb --ensemble2-trajectory filtered_2K39.pdb etc.
for i in $(seq 0 5); do
	cmdline=$cmdline" --ensemble$((i+1))-trajectory ${pdbs[$i]}"
done

# Run encore
echo "Now running: $cmdline"
echo "Results will be written in file: dimred.log"

$PYTHON encore/similarity.py $cmdline &> dimred.log

