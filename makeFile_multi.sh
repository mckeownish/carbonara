#/bin/sh
echo "================================================================================"
echo "Building carbonara..."

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


g++ -c -O3 -std=gnu++14 -o src/main_multi.o src/main_multi.cpp
g++ -O3 -std=gnu++14 -o predictMulti src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/experimentalData.o src/hydrationShellRandom.o src/skmt.o src/writheFP.o src/moleculeFitAndState.o src/main_multi.o src/Logger.o src/helpers.o src/parameters.o

g++ -c -O3 -std=gnu++14 -o src/Flexible_generator.o src/Flexible_generator.cpp

g++ -O3 -std=gnu++14 -o generate_structure src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/Flexible_generator.o

echo "\n"
echo "Built carbonara!"
echo "================================================================================"
echo "Running parallel carbonara..."

echo "\n"
sh RunMe_multi.sh
echo "\n"

echo "Ran parallel carbonara!"
echo "\n"


