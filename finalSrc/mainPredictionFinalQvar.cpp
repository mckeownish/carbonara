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
  argv[ 8] request to apply hydrophobic covering BETWEEN monomers -- currently not used
  argv[ 9] kmin
  argv[10] kmax
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
  
  /* Determine which sections are being altered */
  std::vector< std::vector<int> > vary_sec_list_list;
  determineVaryingSections(argv, vary_sec_list_list);
    
  /* Read in any fixed distances constraints (contact predictions/sulfide bonds) */
  readFixedDistancesConstraints(argv, mol);

  /* Read in the permissible mixture list */
  readPermissibleMixtures(argv, params);
  
  /* Read in the scattering and set up the scattering model */
  experimentalData ed(argv[1]);

  /* Random generator */
  RandomGenerator rng;
 
  /* initialise the state of mol vector */
  moleculeFitAndState molState(mol, params);

  int improvementIndex=0;
  // If we resume from previous run - argv[3] restart True/False
  if((strcmp(argv[3],"True") == 0)){
    improvementIndex=std::atoi(argv[17]);
  }

  std::pair<double,double> overallFit;
  // std::cout << " first getOverallFit \n";
  overallFit = molState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  // if(params.affineTrans==true){
  //     overallFit = molState.getOverallFitForceConnection(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  // } else {
  //     overallFit = molState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  // }
  // std::cout << " first getOverallFit - done \n";

  logger.logMetadata(argv[16], params);

  std::string scatterNameInitial = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, "initial");
  std::string xyzNameInitial = write_molecules(argv[12], improvementIndex, mol, "initial");

  // log starting point
  logger.logEntry(0, 0, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(), 
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameInitial, xyzNameInitial);
  
  logger.consoleInitial(overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(), molState.getDistanceConstraints());
  

  /* Main algorithm */

  /* 

     molState - moleculeFitAndState
     overallFit - molState.getOverallFit() {.first = sum of scat/overlap/dist/writhe, .second = scat only}
     
   */

  // noSections vector tells us how many subsections are in each molecule
  // e.g. for a monomer/dimer mixture noSections[0]=1,noSections[1]=2.
  /* find number of sections in each molecule in mol (vector <ktlMolecule>) */
  std::vector<int> noSections = findNumberSections(mol);

  /* initialise the set of historical states - currently basic, but used to save previous fit stages */
  std::vector<moleculeFitAndState> molStateSet = makeHistoricalStateSet(molState, params);
  
  // loop number
  int k=0;

  // This is a monster while loop - strap in chaps
  while(k < params.noScatterFitSteps){
    
    // Increasing the kmax if we have a good enough fit, consider a little more of the experimental data!
    if(overallFit.second <0.0005 || (params.improvementIndexTest>std::round(params.noScatterFitSteps/5)&& overallFit.second <0.0007)){
      increaseKmax(overallFit, molStateSet, ed, params, logger);
    }
    
    params.improvementIndexTest = params.improvementIndexTest + 1;

    // pick a 'random' molState from the historical molStateSet
    // to become update function
    int index = rng.getChangeIndexProbability(k, params);
    molState = molStateSet[index];
    mol = molState.getMolecule();
    overallFit = molState.getFit();

    /* START TEST PRINT */
    std::cout << " * " << k << " * \n";
    std::cout << "overallFit (updated with improved fitting): " << overallFit.first << "\n";

    moleculeFitAndState TESTmolState = molStateSet[0];
    std::pair<double,double> TESTcurrentFit = TESTmolState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);

    std::cout << "molState from molStateSet - getFit(): " << TESTmolState.getFit().first << "\n";
    std::cout << "molState from molStateSet - getOverallFit(): " << TESTcurrentFit.first << "\n \n";

    std::vector<ktlMolecule> TESTmol = molStateSet[0].getMolecule();
    moleculeFitAndState TESTmolState2(TESTmol, params);
    std::pair<double,double> TESTcurrentFit2 = TESTmolState2.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);

    std::cout << "molState from molStateSet 2 - getFit(): " << TESTmolState2.getFit().first << "\n";
    std::cout << "molState from molStateSet 2 - getOverallFit(): " << TESTcurrentFit2.first << "\n \n";
    /* END TEST PRINT */

    for(int l=0;l<mol.size();l++){

      int netIndex=0;
         
      //loop over the sections of the given molecule (i.e. if its a monomer this loop is tivial, but not for a multimer
      // another monster looooooop
      for(int i=1;i<=noSections[l];i++){
	   
  // Selected transformation option? 
	if(params.affineTrans==true){
        
	  ktlMolecule molCopyR = mol[l];
        
	  double angle = rng.getRotAng(); double theta = rng.getTheAng(); double phi = rng.getPhiAng();
    point kv(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta));

    double xtran = rng.getDistTran(); double ytran = rng.getDistTran(); double ztran = rng.getDistTran();
	  point tranVec(xtran,ytran,ztran);

	  molCopyR.changeMoleculeMultiRotate(angle,kv,i,tranVec);
	  bool cacaDist = molCopyR.checkCalphas(i); 
    
    // Logic here repeated - function!
	  if(cacaDist==false){
	    
	    // calculate the new fit for this
	    moleculeFitAndState newMolState = molState;
	    std::pair<double,double> newOverallFit = newMolState.getOverallFit(ed,params.mixtureList,params.helRatList,molCopyR,params.kmin,params.kmaxCurr,l);
	    double uProb = rng.getDistributionR();
          
	    if(checkTransition(newOverallFit.first, overallFit.first, uProb, k, params.noScatterFitSteps)){

	      improvementIndex++;
        updateAndLog(improvementIndex, mol, molCopyR, molState, newMolState, overallFit, newOverallFit, logger, l, k, ed, params);

        logger.consoleChange("fitImprove", params);
	    } 
	  }
	} // rotate/translate section ends

	// net index tells us how far we are through the whole moelcule
	if(i>1){
	  netIndex=netIndex+mol[l].getSubsecSize(i-1);
	}
	
  // not sure - Chris - all varying sections up for grabs?
  bool doAll = false;

  // Now loop over the secondary structures of the given unit or section
	for(int j=0;j<mol[l].getSubsecSize(i)-1;j++){
	  //std::cout<<" mol "<<l<<" sec "<<i<<" has this many sections "<<mol[l].getSubsecSize(i)<<"\n";
	  int totalIndex = netIndex+j;
        
	  // in this if statement we check which secondary sections are being changed 
	  if((doAll==true) || (std::find(vary_sec_list_list[l].begin(),vary_sec_list_list[l].end(),totalIndex)!=vary_sec_list_list[l].end())){
	    
      // print statement currently in to check what we are changing is correct
	    //std::cout<<" section "<<totalIndex<<" of unit "<<i<<" "<<" sub set number "<<totalIndex-netIndex<<" being altered "<<mol.getSubsecSize(i)<<"\n";
	    
	    // tl;dr - copy the molecule to change it and test if we do better

      int indexCh = totalIndex-netIndex;
      ktlMolecule newMol = mol[l];
      bool cacaDist = modifyMolecule(newMol, mol[l], indexCh, i);

      // Logic here repeated - function!
	    if(cacaDist==false){

	      moleculeFitAndState newmolState = molState;

	      //calculate the fitting of changed molecule
	      std::pair<double,double> newOverallFit = newmolState.getOverallFit(ed,params.mixtureList,params.helRatList,newMol,params.kmin,params.kmaxCurr,l);
	      double uProb = rng.getDistributionR();

	      if(checkTransition(newOverallFit.first, overallFit.first, uProb, k, params.noScatterFitSteps)){

           // Success! Add to the update index
          improvementIndex++;
          updateAndLog(improvementIndex, mol, newMol, molState, newmolState, overallFit, newOverallFit, logger, l, k, ed, params);
          logger.consoleChange("fitImprove", params);
        }
        
	    } 
	  } // if doAll ..
	} // end of j for loop - number of subsections (getSubsecSize)
      } // end of i for loop - noSections in each l
    } // end of l for loop - mol.size()


    // Assign the new 'improved' molecule state to the historical tracker
    molStateSet[index] = molState;
    molStateSet[index].updateMolecule(mol);
    sortVec(molStateSet);
        
    // Print out to terminal window
    logger.consoleFitAttempt(k, improvementIndex, params, overallFit.first, overallFit.second);

    k++;
  } // end of while loop: k < noScatterFitSteps
  
  improvementIndex++;

  // pull the 'best' fit from the historical tracked {remember sorted - 0 index best fitting}
  std::vector<ktlMolecule> molBest = molStateSet[0].getMolecule();

  std::string moleculeNameEnd = write_molecules(argv[12], improvementIndex, mol, "end");
  
  // regenrate molecule hydration layer to update the fit
  moleculeFitAndState molStateBest(molBest, params);

  std::pair<double,double> overallFitBest;
  overallFitBest = molStateBest.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  // if(params.affineTrans==true){
  //     overallFitBest = molStateBest.getOverallFitForceConnection(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  // }else{
  //     overallFitBest= molStateBest.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  // }
  

  // molFitOut.writeScatteringToFile(ed,kmin,kmaxCurr,argv[13]);
  std::string scatterNameEnd = write_scatter(argv[12], improvementIndex, molStateBest, ed, params.kmin, params.kmaxCurr, "end");
  
  std::cout << "\n best overall mol name: " << moleculeNameEnd << "\n";
  std::cout << " overallFitBest fit: " << overallFitBest.first << "\n";

  logger.logEntry(improvementIndex, k, overallFitBest.first, molStateBest.getWrithePenalty(), molStateBest.getOverlapPenalty(), 
                  molStateBest.getDistanceConstraints(), params.kmaxCurr, scatterNameEnd, moleculeNameEnd);

} // end of main