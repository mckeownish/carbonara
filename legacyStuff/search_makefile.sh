#/bin/sh

if ls src/.o 1> /dev/null 2>&1; then
rm src/.o

fi

g++ -c -O3 -std=gnu++14 -o src/point.o  src/point.cpp
g++ -c -O3 -std=gnu++14 -o src/polyHelix.o src/polyHelix.cpp
g++ -c -O3 -std=gnu++14 -o src/randomMolGen.o src/randomMolGen.cpp
g++ -c -O3 -std=gnu++14 -o src/ktlMoleculeRandom.o src/ktlMoleculeRandom.cpp
g++ -c -O3 -std=gnu++14 -o src/experimentalData.o src/experimentalData.cpp
g++ -c -O3 -std=gnu++14 -o src/hydrationShellRandom.o src/hydrationShellRandom.cpp
g++ -c -O3 -std=gnu++14 -o src/skmt.o src/skmt.cpp
g++ -c -O3 -std=gnu++14 -o src/writheFP.o src/writheFP.cpp
g++ -c -O3 -std=gnu++14 -o src/moleculeFitAndState.o src/moleculeFitAndState.cpp
g++ -c -O3 -std=gnu++14 -o src/Logger.o src/Logger.cpp
g++ -c -O3 -std=gnu++14 -o src/helpers.o src/helpers.cpp
g++ -c -O3 -std=gnu++14 -o src/parameters.o src/parameters.cpp


g++ -c -O3 -std=gnu++14 -o src/searchAlgorithm.o src/searchAlgorithm.cpp
g++ -O3 -std=gnu++14 -o searchAlgo src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/experimentalData.o src/hydrationShellRandom.o src/skmt.o src/writheFP.o src/moleculeFitAndState.o src/searchAlgorithm.o src/Logger.o src/helpers.o src/parameters.o