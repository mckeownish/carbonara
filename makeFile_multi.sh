#/bin/sh
echo "================================================================================"
echo "Building carbonara..."

if ls finalSrc/.o 1> /dev/null 2>&1; then
rm finalSrc/.o
fi

g++ -c -O3 -std=gnu++14 -o finalSrc/point.o  finalSrc/point.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/polyHelix.o finalSrc/polyHelix.cpp 
g++ -c -O3 -std=gnu++14 -o finalSrc/randomMolGen.o finalSrc/randomMolGen.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/ktlMoleculeRandom.o finalSrc/ktlMoleculeRandom.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/experimentalData.o finalSrc/experimentalData.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/hydrationShellRandom.o finalSrc/hydrationShellRandom.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/skmt.o finalSrc/skmt.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/writheFP.o finalSrc/writheFP.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/moleculeFitAndState.o finalSrc/moleculeFitAndState.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/Logger.o finalSrc/Logger.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/helpers.o finalSrc/helpers.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/parameters.o finalSrc/parameters.cpp


g++ -c -O3 -std=gnu++14 -o finalSrc/main_multi.o finalSrc/main_multi.cpp
g++ -O3 -std=gnu++14 -o predictMulti finalSrc/point.o finalSrc/polyHelix.o finalSrc/randomMolGen.o finalSrc/ktlMoleculeRandom.o finalSrc/experimentalData.o finalSrc/hydrationShellRandom.o finalSrc/skmt.o finalSrc/writheFP.o finalSrc/moleculeFitAndState.o finalSrc/main_multi.o finalSrc/Logger.o finalSrc/helpers.o finalSrc/parameters.o

g++ -c -O3 -std=gnu++14 -o finalSrc/Flexible_generator.o finalSrc/Flexible_generator.cpp

g++ -O3 -std=gnu++14 -o generate_structure finalSrc/point.o finalSrc/polyHelix.o finalSrc/randomMolGen.o finalSrc/ktlMoleculeRandom.o finalSrc/Flexible_generator.o

echo "\n"
echo "Built carbonara!"
echo "================================================================================"
echo "Running parallel carbonara..."

echo "\n"
sh RunMe_multi.sh
echo "\n"

echo "Ran parallel carbonara!"
echo "\n"


