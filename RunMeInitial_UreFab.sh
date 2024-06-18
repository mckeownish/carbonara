#!/bin/bash
 ScatterFile=/home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/Saxs.dat
 fileLocs=/home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/
 initialCoordsFile=frompdb
 noStructures=1
 pairedPredictions=False
 fixedsections=/home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/varyingSectionSecondary1.dat
 withinMonomerHydroCover=none
 kmin=0.01
 kmax=0.5
 kstart=0.15
 maxNoFitSteps=1
 affineTrans=False
for i in {1..1}

do

   echo " Run number : $i "

   ./getInitialPrediction $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $kmin $kmax $kstart $maxNoFitSteps /home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/tmp/mol$i /home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/tmp/scatter$i.dat /home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/mixtureFile.dat /home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/redundant /home/rdata/ccck67/carbonara-outputs/UreFab/newFitData/tmp/fitLog$i.dat null $affineTrans

done