. ../set_environment.sh
ln -s $ENCORE_BUILD .
ln -s $MDANALYSIS_BUILD . 

# List of files containing the ensembles to be compared. The x-ray derived ensemble is included as a .DCD trajectory file.
ffs=(filtered_1XQQ.pdb  filtered_2K39.pdb  filtered_2KOX.pdb filtered_2LJ5.pdb filtered_2NR2.pdb xray.dcd)

# number of ensembles to be used
nensembles=6

# command line
cmdline="./similarity.py -v --topology filtered_reference.pdb --mode=harmonic --nensembles=$nensembles"

# for each ensemble: add the --ensemble-trajectory option to the command line. For instance: --ensemble1-trajectory filtered_1XQQ.pdb --ensemble2-trajectory filtered_2K39.pdb etc.
for i in $(seq 0 5); do
	cmdline=$cmdline" --ensemble$((i+1))-trajectory ${ffs[$i]}"
done

# Run encore
echo "Now running: $cmdline"
echo "Results will be written in file: harmonic.log"
echo  $ENCORE_BUILD/similarity.py
$PYTHON encore/similarity.py $cmdline
