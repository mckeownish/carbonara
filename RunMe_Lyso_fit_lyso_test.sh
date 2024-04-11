 #!/bin/bash
 ScatterFile=/Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/Saxs.dat
 fileLocs=/Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/
 initialCoordsFile=frompdb
 noStructures=1
 pairedPredictions=False
 fixedsections=/Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/varyingSectionSecondary1.dat
 withinMonomerHydroCover=none
 betweenMonomerHydroCover=none
 kmin=0.01
 kmax=0.25
 maxNoFitSteps=50000
 affineTrans=False
for i in {1..1}

do

   echo " Run number : $i "

   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps /Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/fit_lyso_test/mol$i /Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/fit_lyso_test/scatter$i.dat /Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/mixtureFile.dat /Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/redundant /Users/josh/Documents/PhD/DevDungeon/Feb2024_carbonara/newFitData/Lyso/fit_lyso_test/fitLog$i.dat null $affineTrans

done
