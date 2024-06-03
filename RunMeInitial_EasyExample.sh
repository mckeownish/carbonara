#!/bin/bash
 ScatterFile=newFitData/EasyExample/Saxs.dat
 fileLocs=newFitData/EasyExample/
 initialCoordsFile=frompdb
 noStructures=1
 pairedPredictions=False
 fixedsections=newFitData/EasyExample/varyingSectionSecondary1.dat
 withinMonomerHydroCover=none
 betweenMonomerHydroCover=none
 kmin=0.0
 kmax=1
 maxNoFitSteps=0.25
 affineTrans=False
for i in {1..1}

do

   echo " Run number : $i "

   ./getInitialPrediction $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps newFitData/EasyExample/tmp/mol$i newFitData/EasyExample/tmp/scatter$i.dat newFitData/EasyExample/mixtureFile.dat newFitData/EasyExample/redundant newFitData/EasyExample/tmp/fitLog$i.dat null $affineTrans

done