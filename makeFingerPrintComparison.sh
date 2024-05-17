#/bin/sh

g++ -c -O3 -std=gnu++14 -o src/point.o  src/point.cpp
g++ -c -O3 -std=gnu++14 -o src/polyHelix.o src/polyHelix.cpp
g++ -c -O3 -std=gnu++14 -o src/randomMolGen.o src/randomMolGen.cpp
g++ -c -O3 -std=gnu++14 -o src/ktlMoleculeRandom.o src/ktlMoleculeRandom.cpp
g++ -c -O3 -std=gnu++14 -o src/skmt.o src/skmt.cpp
g++ -c -O3 -std=gnu++14 -o src/writheFP.o src/writheFP.cpp
g++ -c -O3 -std=gnu++14 -o src/mainFileCompareFingerprints.o src/mainFileCompareFingerprints.cpp

g++ -O3 -std=gnu++14 -o compareFingerprints src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/skmt.o src/writheFP.o src/mainFileCompareFingerprints.o
