#!/bin/bash

# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <runX>"
    exit 1
fi

# Get the run folder name from the argument
RUN_FOLDER="$1"

# Determine the root directory based on the script location
ROOT=$(dirname "$(readlink -f "$0")")
# Directory to clear before running
CLEAR_DIR="$ROOT/newFitData/calmodulin_WAXSiS/$RUN_FOLDER"
# Clear the directory
echo "Clearing directory: $CLEAR_DIR"
rm -rf "$CLEAR_DIR"
mkdir -p "$CLEAR_DIR"

### argv[ 1] scattering data file
ScatterFile=$ROOT/newFitData/calmodulin_WAXSiS/Saxs.dat
### argv[ 2] sequence file location
fileLocs=$ROOT/newFitData/calmodulin_WAXSiS/
### argv[ 3] restart tag (use to start from existing prediction)
initialCoordsFile=frompdb
### argv[ 4] paired distances file (can be empty)
pairedPredictions=$ROOT/newFitData/calmodulin_WAXSiS/fixedDistanceConstraints1.dat
echo "\n"
echo "with dist constraints"
echo "\n"
### argv[ 5] fixed sections file (again can be empty)
fixedsections=$ROOT/newFitData/calmodulin_WAXSiS/varyingSectionSecondary1.dat
### argv[ 6] number of structures
noStructures=1
### argv[ 7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it -- Currently not used
withinMonomerHydroCover=none
### argv[ 8] request to apply hydrophobic covering BETWEEN monomers will be a list of pairs to try to hydropobically pair -- currently not used
betweenMonomerHydroCover=none
### argv[ 9] kmin
kmin=0.00
### argv[10] kmax
kmax=0.45
### argv[11] Max number of fitting steps
maxNoFitSteps=25000
### argv[12] prediction file - mol[i] in the fitting folder
predictionFile=$ROOT/newFitData/calmodulin_WAXSiS/$RUN_FOLDER
### argv[13] scattering output file
scatterOut=$ROOT/newFitData/calmodulin_WAXSiS/$RUN_FOLDER
### argv[14] mixture list file, alist of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
mixtureFile=$ROOT/newFitData/calmodulin_WAXSiS/mixtureFile.dat
### argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
prevFitStr=$ROOT/newFitData/calmodulin_WAXSiS/redundant
### argv[16] log file location
logLoc=$ROOT/newFitData/calmodulin_WAXSiS/$RUN_FOLDER
### argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
endLinePrevLog=null
### argv[18] is true if we want to apply affine rotations,false if not.
affineTrans=False

for i in {1..1}
do
    echo "\n"
    echo " >> Run number : $i "
    echo "\n"
    $ROOT/build/bin/predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/fitLog$i.dat $endLinePrevLog $affineTrans
done