#!/bin/bash
 ScatterFile=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/Saxs.dat
 fileLocs=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/
 initialCoordsFile=frompdb
 noStructures=1
 pairedPredictions=False
 fixedsections=/home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/varyingSectionSecondary1.dat
 withinMonomerHydroCover=none
 kmin=0.0
 kmax=0.6
 kstart=0.15
 maxNoFitSteps=1
 affineTrans=False
for i in {1..1}

do

   echo " Run number : $i "

   ./getInitialPrediction $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $kmin $kmax $kstart $maxNoFitSteps /home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/tmp/mol$i /home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/tmp/scatter$i.dat /home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/mixtureFile.dat /home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/redundant /home/b21user/Documents/carbonara/newFitData/Lysozyme/newFitData/tmp/fitLog$i.dat null $affineTrans

done