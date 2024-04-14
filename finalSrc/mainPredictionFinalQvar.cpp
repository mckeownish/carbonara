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

// todo: 
// add helRatList, mixtureList to params (i think)
// .getOverallFit (less args)
// improved cout - logger function, much easier to read on terminal
// test with cameron multi-mol - work? eeesh
 // rotate/translate become own function?
 // .logEntry (params / less)

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

  /* Read in the molecule(s) */
  std::vector<ktlMolecule> mol;
  readInStructures(argv, mol, params);
  
  /* Determine which sections are being altered */
  std::vector< std::vector<int> > vary_sec_list_list;
  determineVaryingSections(argv, vary_sec_list_list);
    
  /* Read in any fixed distances constraints (contact predictions/sulfide bonds) */
  readFixedDistancesConstraints(argv, mol);

  /* Read in the permissible mixture list */
  std::vector< std::vector<double> > mixtureList;
  readPermissibleMixtures(argv, mixtureList);
  
  /* Read in the scattering and set up the scattering model */
  experimentalData ed(argv[1]);

  /* Random generator */
  RandomGenerator rng;
 
  /* initialise the first fit */
  moleculeFitAndState molFit(mol, params);

  std::pair<double,double> scatterFit;

  int improvementIndex=0;
  // If we resume from previous run - argv[3] restart True/False
  if((strcmp(argv[3],"True") == 0)){
    improvementIndex=std::atoi(argv[17]);
  }
  
  // redo molecule hydration layer to update fit
  // ** moleculeFitAndState **
  moleculeFitAndState molFitOg(mol, params);
  std::pair<double,double> scatterFitOut = molFitOg.getOverallFit(ed, mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  std::string scatterNameInitial = write_scatter(argv[12], improvementIndex, molFitOg, ed, params.kmin, params.kmaxCurr, "initial");

  auto start = high_resolution_clock::now();
  auto curr = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(curr - start);

  logger.logMetadata(argv[16], params, params.helRatList);
  
  std::string tempMolName = "temp_molecule_name";
  // log starting point
  logger.logEntry(0, 0, scatterFitOut.first, molFitOg.getWrithePenalty(), molFitOg.getOverlapPenalty(), 
                  molFitOg.getDistanceConstraints(), duration.count(), params.kmaxCurr, scatterNameInitial, tempMolName);
  
  logger.consoleInitial(scatterFitOut.first, molFitOg.getWrithePenalty(), molFitOg.getOverlapPenalty(), molFitOg.getDistanceConstraints());
    
  /* Main algorithm */
  // This whole section is bracket hell
  // My god - we need break stuff tf up into smaller functions
  
  // noSections vector tells us how many subsections are in each molecule
  // e.g. for a monomer/dimer mixture noSections[0]=1,noSections[1]=2.
  std::vector<int> noSections;
    
  for(int i=0;i<mol.size();i++){
    int noSectionsTmp = mol[i].noChains();
    noSections.push_back(noSectionsTmp);
  }
    
  // argv[11] Max number of fitting steps
  int noScatterFitSteps=std::atoi(argv[11]);

  // set up the vector of existing structures
  // initialise to the initial guess

  int noHistoricalFits = 1;
  std::vector<moleculeFitAndState> molFitAndStateSet;
  
  for(int i=0;i<noHistoricalFits;i++){
    molFitAndStateSet.push_back(molFit);
  }

  int improvementIndexTest=0;
  
  // loop number
  int k=0;

  // This is a monster while loop - strap in chaps
  while(k < noScatterFitSteps){
      
    // Print out to terminal window
    logger.consoleFitAttempt(k, improvementIndex, params, scatterFit.first, scatterFit.second);
    // std::cout<<params.kmin<<" "<<params.kmaxCurr<<" "<<params.kmax<<" "<<scatterFit.first<<" "<<scatterFit.second<<"\n";
    
    // Increasing the kmax if we have a good enough fit, consider a little more of the experimental data!
    if(scatterFit.second <0.0005 || (improvementIndexTest>std::round(noScatterFitSteps/5)&& scatterFit.second <0.0007)){
      
      // if we have achieved a sufficiently good fit include more data.
      params.kmaxCurr=params.kmaxCurr+0.01;
        
      if(params.kmaxCurr>params.kmax){
        params.kmaxCurr=params.kmax;
      }
      
      logger.consoleChange("krangeIncrease", params);
        
      improvementIndexTest=0;
      // generate a new first fit.
      scatterFit = molFitAndStateSet[0].getOverallFit(ed,mixtureList,params.helRatList,params.kmin,params.kmaxCurr);
    }
    
    improvementIndexTest++;
    // loop over the molecules (e.g monomer and dimer fit

    int index = rng.getChangeIndexProbability(k, noHistoricalFits, noScatterFitSteps);
    molFit = molFitAndStateSet[index];
    mol = molFit.getMolecule();
    scatterFit = molFit.getFit();

    for(int l=0;l<mol.size();l++){

      int netIndex=0;
         
      //loop over the sections of the given molecule (i.e. if its a monomer this loop is tivial, but not for a multimer
      
      // another monster looooooop
      for(int i=1;i<=noSections[l];i++){
	
    // To become own function?
    // rotate/translate section begins
	if(params.affineTrans==true){
        
    // bool transformationApplied = performAffineTransformation(mol, l, generator1, distTran, rotAng, theAng, phiAng, 
    //                                                      distributionR, generator, molFit, ed, mixtureList, 
    //                                                      helRatList, kmin, kmaxCurr, k, noScatterFitSteps);

	  ktlMolecule molCopyR = mol[l];
        
	  double angle = rng.getRotAng(); double theta = rng.getTheAng(); double phi = rng.getPhiAng();
    point kv(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta));

    double xtran = rng.getDistTran(); double ytran = rng.getDistTran(); double ztran = rng.getDistTran();
	  point tranVec(xtran,ytran,ztran);

	  molCopyR.changeMoleculeMultiRotate(angle,kv,i,tranVec);
	  bool cacaDist= molCopyR.checkCalphas(i); 
        
	  if(cacaDist==false){
	    
	    // calculate the new fit for this
	    moleculeFitAndState molFitTmp = molFit ;
	    //calculate all amino acid distances for changed molecule
	    std::pair<double,double> fitTemp = molFitTmp.getOverallFit(ed,mixtureList,params.helRatList,molCopyR,params.kmin,params.kmaxCurr,l);
	    // check if we have imporved
	    //std::cout<<"how change ? "<<fitTemp<<" "<<scatterFit<<"\n";
	    double uProb = rng.getDistributionR();;
          
	    if(checkTransition(fitTemp.first,scatterFit.first,uProb,k,noScatterFitSteps)){

	      scatterFit = fitTemp;
	      mol[l] = molCopyR;
	      molFit = molFitTmp;

	      // to output during fitting to "show the process"
	      improvementIndex++;

        std::string moleculeNameTrans = write_molecules(argv[12], improvementIndex, mol);
        std::string scatterNameTrans = write_scatter(argv[12], improvementIndex, molFitTmp, ed, params.kmin, params.kmaxCurr);

	      curr = high_resolution_clock::now();
	      duration = duration_cast<microseconds>(curr - start);

        // log file write
        logger.logEntry(improvementIndex, k, scatterFit.first, molFitTmp.getWrithePenalty(), molFitTmp.getOverlapPenalty(), 
                        molFitTmp.getDistanceConstraints(), duration.count(), params.kmaxCurr, scatterNameTrans, moleculeNameTrans);

        logger.consoleChange("fitImprove", params);

	    } // transition check end
	  } // cacaDist end§
	} // rotate/translate section ends

	
	// net index tells us how far we are through the whole moelcule
	if(i>1){
	  netIndex=netIndex+mol[l].getSubsecSize(i-1);
	}
	
	// Now loop over the secondary structures of the given unit or section

  // not sure - Chris - all varying sections up for grabs?
  bool doAll = false;


	for(int j=0;j<mol[l].getSubsecSize(i)-1;j++){
	  //std::cout<<" mol "<<l<<" sec "<<i<<" has this many sections "<<mol[l].getSubsecSize(i)<<"\n";
	  int totalIndex = netIndex+j;
        
	  // in this if statement we check which secondary sections are being changed 
	  if((doAll==true) || (std::find(vary_sec_list_list[l].begin(),vary_sec_list_list[l].end(),totalIndex)!=vary_sec_list_list[l].end())){
	    // print statement currently in to check what we are changing is correct
	    //std::cout<<" section "<<totalIndex<<" of unit "<<i<<" "<<" sub set number "<<totalIndex-netIndex<<" being altered "<<mol.getSubsecSize(i)<<"\n";
	    
	    // tl;dr - copy the molecule to change it and test if we do better

	    ktlMolecule molCopy = mol[l];

	    int indexCh = totalIndex-netIndex;
	    molCopy.changeMoleculeSingleMulti(indexCh,i);

	    // this (checkCalphas) checks if there haven't been any rouge sections created (some occasional flaws in the procedure which are to be ironed out
	    bool cacaDist= molCopy.checkCalphas(i,mol[l]);
	    if(cacaDist==false){

	      // calculate the new fit for this
	      moleculeFitAndState molFitTmp = molFit;

	      //calculate all amino acid distances for changed molecule
	      std::pair<double,double> fitTemp = molFitTmp.getOverallFit(ed,mixtureList,params.helRatList,molCopy,params.kmin,params.kmaxCurr,l);

	      // check if we have imporved
	      //std::cout<<"how change ? "<<fitTemp<<" "<<scatterFit<<"\n";
	      double uProb = rng.getDistributionR();

	      if(checkTransition(fitTemp.first,scatterFit.first,uProb,k,noScatterFitSteps)){

          // the ol' update shuffle
          scatterFit = fitTemp;
          mol[l] = molCopy;
          molFit = molFitTmp;

          // Success! Add to the update index
          improvementIndex++;

          // writing functs from helpers.cpp
          std::string moleculeNameMain = write_molecules(argv[12], improvementIndex, mol);
          std::string scatterNameMain = write_scatter(argv[12], improvementIndex, molFitTmp, ed, params.kmin, params.kmaxCurr);

          curr = high_resolution_clock::now();
          duration = duration_cast<microseconds>(curr - start);

          logger.logEntry(improvementIndex, k, scatterFit.first, molFitTmp.getWrithePenalty(), molFitTmp.getOverlapPenalty(), 
                          molFitTmp.getDistanceConstraints(), duration.count(), params.kmaxCurr, scatterNameMain, moleculeNameMain);

	      } // check transition end
            
	    } // if cacaDist == False (after making a random change) end
	  } // if doAll &

	} // end of j for loop - number of subsections (getSubsecSize)

      } // end of i for loop - noSections in each l
    } // end of l for loop - mol.size()

    molFitAndStateSet[index] = molFit;
    molFitAndStateSet[index].updateMolecule(mol);
    sortVec(molFitAndStateSet);
      
    for(int i=0;i<noHistoricalFits;i++){
      // std::cout<<"step "<<k<<" "<<i<<" "<<molFitAndStateSet[i].currFit<<"\n";
      // logger.consoleCurrentStep(k, i, molFitAndStateSet[i].currFit);
    }
    
    // increase k - eventually break the while loop!
    k++;

  } // end of while loop: k < noScatterFitSteps
  
  improvementIndex++;

  mol = molFitAndStateSet[0].getMolecule();

  std::string moleculeNameEnd = write_molecules(argv[12], improvementIndex, mol);
    
  // regenrate molecule hydration layer to update the fit
  moleculeFitAndState molFitOut(mol, params);
  scatterFitOut = molFitOut.getOverallFit(ed, mixtureList, params.helRatList, params.kmin, params.kmaxCurr);

  // molFitOut.writeScatteringToFile(ed,kmin,kmaxCurr,argv[13]);
  std::string scatterNameEnd = write_scatter(argv[12], improvementIndex, molFitOut, ed, params.kmin, params.kmaxCurr, "end");

  logger.logEntry(improvementIndex, k, scatterFit.first, molFitOut.getWrithePenalty(), molFitOut.getOverlapPenalty(), 
                  molFitOut.getDistanceConstraints(), duration.count(), params.kmaxCurr, scatterNameEnd, moleculeNameEnd);

} // end of main
      
