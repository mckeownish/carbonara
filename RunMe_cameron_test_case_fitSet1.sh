#!/bin/bash
 ScatterFile=/Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/Saxs.dat
 fileLocs=/Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/
 initialCoordsFile=frompdb
 noStructures=1
 pairedPredictions=True
 fixedsections=/Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/varyingSectionSecondary1.dat
 withinMonomerHydroCover=none
 betweenMonomerHydroCover=none
 kmin=0.01
 kmax=0.25
 maxNoFitSteps=5000
 affineTrans=True
for i in {1..1}

do

   echo " Run number : $i "

   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps /Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/fitSet1/mol$i /Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/fitSet1/scatter$i.dat /Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/mixtureFile.dat /Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/redundant /Users/josh/Documents/PhD/DevDungeon/March2024_carbonara/newFitData/cameron_test_case/fitSet1/fitLog$i.dat null $affineTrans

done
