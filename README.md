
# CarbonaraPy

CarbonaraPy is a Python package that interfaces with our C++ package Carbonara that provides tools for correcting structural predictions of proteins using X-ray small-angle scattering (SAXS) in solution.


## Features

- Correct structural predictions of proteins based on SAXS data.
- Refinement of crystallograhy or AlphaFold predictions (even low confidence structures!)
- Model mixtures of different structural states, help to identify dynamics in solution
- Easy-to-use Python wrapper.


## Installation

Ensure you have g++ installed, as the package requires compilation of C++ source code.

```bash
pip install carbonarapy
```


## Usage

Here's a simple example of how to use CarbonaraPy:

```python
from carbonarapy as run_predict_structure

# Example function call
run_predict_structure( ScatterFile, fileLocs, initialCoordsFile, pairedPredictions, fixedsections, noStructures,
					   withinMonomerHydroCover, betweenMonomerHydroCover, kmin, kmax, maxNoFitSteps, predictionFile,
					   scatterOut, mixtureFile, prevFitStr, logLoc, endLinePrevLog, affineTrans, 
					   number_runs=3, verbose=True
					 )

```

See the tutorial notebook for help with function arguments! 


## Contributing

Contributions to CarbonaraPy are welcome! Please read our contributing guidelines for more information on how to report issues, submit changes, and contribute to the development.

License

CarbonaraPy is licensed under the MIT License. See the LICENSE file for more details.


## Acknowledgments

This project is developed and maintained by Joshua McKeown, Arron Bale & Christopher Prior @ Durham University. We appreciate all the contributions from the community! 