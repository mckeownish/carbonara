This codebase is a C++ project focused on refinement of protein structures against experimental in solution SAXS. It includes functionality for generating and manipulating molecular structures, analyzing their properties, and fitting them to experimental data.

## Building with CMake

To build the project using CMake, follow these steps:


1. Make sure you have CMake installed on your system (version 3.10 or higher is recommended)

```cmake -version```

2. Open a terminal and navigate to the carbonara root directory:

```cd path/to/carbonara```

3. Create a build directory and navigate into it:

```
mkdir build
cd build
```

4. Generate the build files:

```
cmake ..
```

5. Build the project:

```
make
```

## Key Components

### 1. ktlMolecule Class

This class represents a molecular structure. It includes methods for:
- Reading in sequence and coordinate data
- Manipulating the molecular structure
- Analyzing properties like hydrophobicity and coiled-coil potential

Key methods to focus on:
- `readInSequence()`: Parses sequence data
- `readInCoordinates()`: Loads coordinate data
- `changeMoleculeSingleMulti()`: Modifies a specific part of the molecule

### 2. hydrationShellMinimal Class

Handles the calculation of hydration shells around molecules. Key areas:
- Generation of hydration layer
- Calculation of solvent-molecule distances

### 3. experimentalData Class

Deals with experimental scattering data and fitting. Important methods:
- `fitToScattering()`: Fits molecular model to scattering data
- `setPhases()`: Sets up scattering phases for calculations

### 4. randomMol Class

Generates random molecular structures. Key functionality:
- Creation of random sections with specific properties
- Blending different structural elements (e.g., loops to helices)

### 5. writheFP Class

Calculates writhe (a topological property) for molecular structures. 
- `DIDownSample()`: Calculates downsampled writhe
- `compareFingerPrints()`: Compares writhe "fingerprints" between structures

## Main Algorithms

### Structure Generation and Manipulation

Located primarily in `randomMol` class. The main method to focus on is `makeRandomMolecule()`.

### Scattering Data Fitting

Implemented in `experimentalData` class. Key method is `fitToScatteringMultiple()`.

### Writhe Calculation

Implemented in `writheFP` class. The main method is `DIDownSample()`.

## Main Execution Flow

The primary execution flow is in `mainPredictionFinalQvar.cpp`. It follows these steps:

1. Initialize parameters and data structures
2. Load experimental data
3. Generate or load initial molecular structures
4. Iteratively modify structures and evaluate fit
5. Output results

## Areas for Improvement

1. Code Organization: Many functions, especially in main files, are very long and could be broken down.
2. Error Handling: More robust error checking and handling is needed throughout.
3. Memory Management: Consider replacing raw pointers with smart pointers.
4. Parallelism: There's potential for more parallelism in computationally intensive parts.
5. Testing: Implement unit tests for key components.

## Next Steps for Development

1. ~~Refactor `mainPredictionFinalQvar.cpp` to improve readability and maintainability.~~
2. Implement more comprehensive error handling.
3. Optimize performance-critical sections, possibly using parallel computing techniques.
4. Improve documentation throughout the codebase.
5. Implement a testing framework and write unit tests for key components.
