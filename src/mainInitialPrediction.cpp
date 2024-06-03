/* Carbonara Version: 0.1.9 */

#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h>
#include "moleculeFitAndState.h"
#include <cstring>
#include <chrono>
#include <tuple>

#include "Logger.h"
#include "helpers.h"

using namespace std::chrono;

// note: this version showing funky behaviour with getFit()'s - not always consistent!

/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

  argv[ 1] scattering data file
  argv[ 2] sequence file location
  argv[ 3] restart tag (use to start from existing prediction)
  argv[ 4] paired distances file (can be empty)
  argv[ 5] fixed sections file (again can be empty)
  argv[ 6] number of structures
  argv[ 7] request to apply hydrophobic covering WITHIN monomers -- Currently not used
  argv[ 8] kmin
  argv[ 9] kmax
  argv[10] kstart
  argv[11] Max number of fitting steps
  argv[12] prediction file - mol[i] in the fitting folder
  argv[13] scattering output file
  argv[14] mixture list file, a list of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
  argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
  argv[16] log file location
  argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
  argv[18] is true if we want to apply affine rotations, false if not.

 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- */

int main( int argc, const char* argv[] ) {

  /* initialise the log file */
  Logger logger(argv[16]);

  /* Set up model parameters */
  ModelParameters params = loadParameters(argv);

  /* >> determine initial model: Two options no initial prediction, we must generate a structure
   or some initial structure provided. Actually we need a half-half option */

  /* Initialise the molecule(s) vector */
  std::vector<ktlMolecule> mol;


  readInStructures(argv, mol, params);
    
  /* Read in any fixed distances constraints (contact predictions/sulfide bonds) */
  readFixedDistancesConstraints(argv, mol);
  
  /* Read in the permissible mixture list */
  readPermissibleMixtures(argv, params);

  
  /* Read in the scattering and set up the scattering model */
  experimentalData ed(argv[1]);
  
  /* Random generator */
  RandomGenerator rng;
 
  moleculeFitAndState molState(mol, params);

  std::pair<double,double> overallFit;

  int improvementIndex=0;
  
  if(params.affineTrans==true){
      overallFit = molState.getOverallFitForceConnection(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  }else {
      overallFit = molState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  }
  // std::cout << " first getOverallFit - done \n";


  std::string scatterNameInitial = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, "initial");
  std::string xyzNameInitial = write_molecules(argv[12], improvementIndex, mol, "initial");
 
} // end of main
