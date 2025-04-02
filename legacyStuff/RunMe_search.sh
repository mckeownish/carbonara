#!/bin/bash

# Determine the root directory based on the script location
ROOT=$(dirname "$(readlink -f "$0")")

 ### argv[ 1] scattering data file
 ScatterFile=$ROOT/newFitData/test_simple/Saxs.dat

 ### argv[ 2] sequence file location
 fileLocs=$ROOT/newFitData/test_simple/

 ### argv[ 3] restart tag (use to start from existing prediction)
 initialCoordsFile=frompdb

 ### argv[ 4] paired distances file (can be empty)
 pairedPredictions=False

 ### argv[ 5] fixed sections file (again can be empty)
 fixedsections=$ROOT/newFitData/test_simple/varyingSectionSecondary1.dat

 ### argv[ 6] number of structures
 noStructures=1

 ### argv[ 7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it -- Currently not used
 withinMonomerHydroCover=none

 ### argv[ 8] request to apply hydrophobic covering BETWEEN monomers will be a list of pairs to try to hydropobically pair -- currently not used
 betweenMonomerHydroCover=none

 ### argv[ 9] kmin
 kmin=0.01

 ### argv[10] kmax
 kmax=0.25

 ### argv[11] Max number of fitting steps
 maxNoFitSteps=50

 ### argv[12] prediction file - mol[i] in the fitting folder
 predictionFile=$ROOT/newFitData/test_simple/fitdata

 ### argv[13] scattering output file
 scatterOut=$ROOT/newFitData/test_simple/fitdata

 ### argv[14] mixture list file, alist of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
 mixtureFile=$ROOT/newFitData/test_simple/mixtureFile.dat

 ### argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
 prevFitStr=$ROOT/newFitData/test_simple/redundant

 ### argv[16] log file location
 logLoc=$ROOT/newFitData/test_simple/fitdata

 ### argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
 endLinePrevLog=null

 ### argv[18] is true if we want to apply affine rotations,false if not.
 affineTrans=False

for i in {1..1}

do
    echo "\n"
    echo " >> Run number : $i "
    echo "\n"
    
    ./searchAlgo $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/searchLog$i.dat $endLinePrevLog $affineTrans
        

done