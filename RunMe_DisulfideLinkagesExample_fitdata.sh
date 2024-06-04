#!/bin/bash
 ScatterFile=newFitData/DisulfideLinkagesExample/Saxs.dat
 fileLocs=newFitData/DisulfideLinkagesExample/
 initialCoordsFile=frompdb
 pairedPredictions=False
 fixedsections=newFitData/DisulfideLinkagesExample/varyingSectionSecondary1.dat
 noStructures=1
 withinMonomerHydroCover=none
 kmin=0.01
 kmax=0.15
 kstart=0.25
 maxNoFitSteps=10000
 predictionFile=newFitData/DisulfideLinkagesExample/fitdata
 scatterOut=newFitData/DisulfideLinkagesExample/fitdata
 mixtureFile=newFitData/DisulfideLinkagesExample/mixtureFile.dat
 prevFitStr=newFitData/DisulfideLinkagesExample/redundant
 logLoc=newFitData/DisulfideLinkagesExample/fitdata
 endLinePrevLog=null
 affineTrans=False
for i in {1..5}

do

   echo " Run number : $i "

   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $kmin $kmax $kstart $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/fitLog$i.dat $endLinePrevLog $affineTrans

done