This codebase is a C++ project focused on refinement of protein structures against experimental in solution SAXS. It includes functionality for generating and manipulating molecular structures, analyzing their properties, and fitting them to experimental data.

## Building with CMake

To build the project using CMake, follow these steps:


1. Open a terminal and make sure you have CMake installed on your system (version 3.10 or higher is recommended)

```
cmake -version
```

2. Navigate to the carbonara root directory:

```
cd path/to/carbonara
```

3. Inside the carbonara directory, create a build directory and navigate into it:

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

## Testing Carbonara

To ensure your version of Carbonara is running, run the test case:

1. Ensure you are located in `/path/to/carbonara` and run the following command

```
sh RunMe_humanSMARCAL1.sh
```

## Using Carbonara for new structures

To refine protein structure predictions with your own SAXS data, you'll need:

1. A PDB starting model (AlphaFold or crystal structure recommended)
2. SAXS experimental data in Å units with three columns: q, I, and I error

### Setting up the Python environment

```bash
# Create a new conda environment
conda create -n carbonara_py python=3.10
conda activate carbonara_py

# Install required packages
pip install pandas 
pip install numpy 
pip install cython 
pip install tqdm 
pip install mdtraj 
pip install biobox
```

Then run:

```bash
python setup_carbonara.py --pdb path/to/pdb --saxs path/to/saxs --name ProteinName 
```
```bash
# Optional flags for customizing refinement
--fit_n_times INT     Number of times to run the fit (default: 5)
--min_q FLOAT         Minimum q-value (default: 0.01)
--max_q FLOAT         Maximum q-value (default: 0.2)
--max_fit_steps INT   Maximum number of fitting steps (default: 1000)
--pairedQ             Use paired predictions
--rotation            Apply affine rotations

```

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
