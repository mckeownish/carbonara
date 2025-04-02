#!/bin/bash

# Determine the root directory based on the script location
ROOT=$(dirname "$(readlink -f "$0")")

# Build directory
BUILD_DIR="$ROOT/build"

# Ensure the build directory exists and the project is built
if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    cmake ..
    make
    cd "$ROOT"
elif [ ! -f "$BUILD_DIR/bin/predictStructureQvary" ]; then
    cd "$BUILD_DIR"
    make
    cd "$ROOT"
fi

# Set up variables (unchanged from your original script)
ScatterFile=$ROOT/newFitData/test_simple/Saxs.dat
fileLocs=$ROOT/newFitData/test_simple/
initialCoordsFile=frompdb
pairedPredictions=False
fixedsections=$ROOT/newFitData/test_simple/varyingSectionSecondary1.dat
noStructures=1
withinMonomerHydroCover=none
betweenMonomerHydroCover=none
kmin=0.01
kmax=0.25
maxNoFitSteps=50
predictionFile=$ROOT/newFitData/test_simple/fitdata
scatterOut=$ROOT/newFitData/test_simple/fitdata
mixtureFile=$ROOT/newFitData/test_simple/mixtureFile.dat
prevFitStr=$ROOT/newFitData/test_simple/redundant
logLoc=$ROOT/newFitData/test_simple/fitdata
endLinePrevLog=null
affineTrans=False

for i in {1..1}
do
    echo -e "\n"
    echo " >> Run number : $i "
    echo -e "\n"

    "$BUILD_DIR/bin/predictStructureQvary" \
        "$ScatterFile" \
        "$fileLocs" \
        "$initialCoordsFile" \
        "$pairedPredictions" \
        "$fixedsections" \
        "$noStructures" \
        "$withinMonomerHydroCover" \
        "$betweenMonomerHydroCover" \
        "$kmin" \
        "$kmax" \
        "$maxNoFitSteps" \
        "$predictionFile/mol$i" \
        "$scatterOut/scatter$i.dat" \
        "$mixtureFile" \
        "$prevFitStr" \
        "$logLoc/fitLog$i.dat" \
        "$endLinePrevLog" \
        "$affineTrans"
done
