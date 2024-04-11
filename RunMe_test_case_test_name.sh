#!/bin/bash
 ROOT = /Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/

 ### argv[ 1] scattering data file
 ScatterFile=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/Saxs.dat

 ### argv[ 2] sequence file location
 fileLocs=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/

 ### argv[ 3] restart tag (use to start from existing prediction)
 initialCoordsFile=frompdb

 ### argv[ 4] paired distances file (can be empty)
 pairedPredictions=True

 ### argv[ 5] fixed sections file (again can be empty)
 fixedsections=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/varyingSectionSecondary1.dat

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
 maxNoFitSteps=500

 ### argv[12] prediction file - mol[i] in the fitting folder
 predictionFile=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/test_name/mol$i

 ### argv[13] scattering output file
 scatterOut=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/test_name/scatter$i.dat

 ### argv[14] mixture list file, alist of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
 mixtureFile=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/mixtureFile.dat

 ### argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
 prevFitStr=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/redundant

 ### argv[16] log file location
 logLoc=/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/test_case/fitLog$i.dat

 ### argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
 endLinePrevLog=null

 ### argv[18] is true if we want to apply affine rotations,false if not.
 affineTrans=True

for i in {1..2}

do

   echo " Run number : $i "

   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps $predictionFile $scatterOut $mixtureFile $prevFitStr $logLoc $endLinePrevLog $affineTrans

done