#!/bin/bash

# Function to print usage information
print_usage() {
    echo "Usage: $0 <scattering_file> <fingerprint_file> <coordinate_file> <output_prefix>"
    echo "  <scattering_file>  : Path to the scattering data file"
    echo "  <fingerprint_file> : Path to the fingerprint (sequence) file"
    echo "  <coordinate_file>  : Path to the coordinate file"
    echo "  <output_prefix>    : Prefix for output files"
}

# Check if correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Error: Incorrect number of arguments."
    print_usage
    exit 1
fi

# Assign command-line arguments to variables
ScatterFile="$1"
FingerprintFile="$2"
CoordinateFile="$3"
OutputPrefix="$4"

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
elif [ ! -f "$BUILD_DIR/bin/single_fit" ]; then
    cd "$BUILD_DIR"
    make
    cd "$ROOT"
fi

# Check if input files exist
for file in "$ScatterFile" "$FingerprintFile" "$CoordinateFile"; do
    if [ ! -f "$file" ]; then
        echo "Error: File not found: $file"
        exit 1
    fi
done

# Run the single_fit program
"$BUILD_DIR/bin/single_fit" \
    "$ScatterFile" \
    "$FingerprintFile" \
    "$CoordinateFile" \
    "$OutputPrefix"

# Print a message indicating the script has finished
echo "Single fit analysis complete. Check the output files with prefix: $OutputPrefix"
