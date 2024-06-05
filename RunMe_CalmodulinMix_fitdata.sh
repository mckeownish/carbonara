#!/bin/bash
 ScatterFile=newFitData/CalmodulinMix/Saxs.dat
 fileLocs=newFitData/CalmodulinMix/
 initialCoordsFile=frompdb
 pairedPredictions=False
 fixedsections=newFitData/CalmodulinMix/varyingSectionSecondary1.dat
 noStructures=2
 withinMonomerHydroCover=none
 kmin=0.01
 kmax=0.34
 kstart=0.15
 maxNoFitSteps=10000
 predictionFile=newFitData/CalmodulinMix/fitdata
 scatterOut=newFitData/CalmodulinMix/fitdata
 mixtureFile=newFitData/CalmodulinMix/mixtureFile.dat
 prevFitStr=newFitData/CalmodulinMix/redundant
 logLoc=newFitData/CalmodulinMix/fitdata
 endLinePrevLog=null
 affineTrans=False
for i in {1..10}

do

   echo " Run number : $i "

   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $kmin $kmax $kstart $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/fitLog$i.dat $endLinePrevLog $affineTrans

done