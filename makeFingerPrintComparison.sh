#/bin/sh

g++ -c -O3 -std=gnu++14 -o finalSrc/point.o  finalSrc/point.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/polyHelix.o finalSrc/polyHelix.cpp 
g++ -c -O3 -std=gnu++14 -o finalSrc/randomMolGen.o finalSrc/randomMolGen.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/ktlMoleculeRandom.o finalSrc/ktlMoleculeRandom.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/skmt.o finalSrc/skmt.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/writheFP.o finalSrc/writheFP.cpp	
g++ -c -O3 -std=gnu++14 -o finalSrc/mainFileCompareFingerprints.o finalSrc/mainFileCompareFingerprints.cpp

g++ -O3 -std=gnu++14 -o compareFingerprints finalSrc/point.o finalSrc/polyHelix.o finalSrc/randomMolGen.o finalSrc/ktlMoleculeRandom.o finalSrc/skmt.o finalSrc/writheFP.o finalSrc/mainFileCompareFingerprints.o
