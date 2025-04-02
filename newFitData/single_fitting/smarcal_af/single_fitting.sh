#!/bin/bash

SINGLE_FIT_BIN="/Users/josh/Documents/PhD/DevDungeon/carbonara/build/bin/single_fit"


CURRENT_DIR=$(pwd)
OUTPUT_PREFIX="fitdata/singlefit"

SCATTER_FILE="Saxs.dat"
FP_FILE="fingerPrint.dat"
COORD_FILE="coordinates1.dat"


# echo "Current directory: $CURRENT_DIR"
# echo "Single_fit executable: $SINGLE_FIT_BIN"



# Normal run
# echo "Running single_fit normally..."
"$SINGLE_FIT_BIN" "${CURRENT_DIR}/${SCATTER_FILE}" "${CURRENT_DIR}/${FP_FILE}" "${CURRENT_DIR}/${COORD_FILE}" "$OUTPUT_PREFIX" 0.15
