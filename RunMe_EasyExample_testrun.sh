#!/bin/bash
 ScatterFile=newFitData/EasyExample/Saxs.dat
 fileLocs=newFitData/EasyExample/
 initialCoordsFile=frompdb
 pairedPredictions=True
 fixedsections=newFitData/EasyExample/varyingSectionSecondary1.dat
 noStructures=1
 withinMonomerHydroCover=none
 betweenMonomerHydroCover=none
 kmin=0.01
 kmax=0.2
 maxNoFitSteps=20
 predictionFile=newFitData/EasyExample/testrun
 scatterOut=newFitData/EasyExample/testrun
 mixtureFile=newFitData/EasyExample/mixtureFile.dat
 prevFitStr=newFitData/EasyExample/redundant
 logLoc=newFitData/EasyExample/testrun
 endLinePrevLog=null
 affineTrans=False
for i in {1..10}

do

   echo " Run number : $i "

   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/fitLog$i.dat $endLinePrevLog $affineTrans

done