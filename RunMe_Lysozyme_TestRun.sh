#!/bin/bash
 ScatterFile=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/Saxs.dat
 fileLocs=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/
 initialCoordsFile=frompdb
 pairedPredictions=False
 fixedsections=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/varyingSectionSecondary1.dat
 noStructures=1
 withinMonomerHydroCover=none
 kmin=0.0
 kmax=0.6
 kstart=0.15
 maxNoFitSteps=10000
 predictionFile=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/TestRun
 scatterOut=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/TestRun
 mixtureFile=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/mixtureFile.dat
 prevFitStr=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/redundant
 logLoc=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/TestRun
 endLinePrevLog=null
 affineTrans=True
for i in {1..10}

do

   echo " Run number : $i "

   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $kmin $kmax $kstart $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/fitLog$i.dat $endLinePrevLog $affineTrans

done