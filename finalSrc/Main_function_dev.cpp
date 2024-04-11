// std::tuple<bool, ktlMolecule, std::pair<double, double>, moleculeFitAndState>
std::tuple<bool, ktlMolecule, std::pair<double, double>, moleculeFitAndState>
performAffineTransformation(ktlMolecule& mol, 
                            int sectionIndex, 
                            const std::default_random_engine& generator1, 
                            const std::uniform_real_distribution<>& distTran, 
                            const std::uniform_real_distribution<>& rotAng, 
                            const std::uniform_real_distribution<>& theAng, 
                            const std::uniform_real_distribution<>& phiAng, 
                            const std::uniform_real_distribution<double>& distributionR, 
                            const std::default_random_engine& generator, 
                            moleculeFitAndState& molFit, 
                            const std::vector<double>& ed, 
                            const std::vector<double>& mixtureList, 
                            const std::vector<double>& helRatList, 
                            double kmin, 
                            double kmaxCurr, 
                            int k, 
                            int noScatterFitSteps) {

  std::pair<double, double> scatterFit; // This should be initialized to the current scatter fit value before calling the function
  
  ktlMolecule molCopyR = mol;

  // all the ingredients for making a random rigid body transformation to our molecule!
  double angle = rotAng(generator1); 
  double theta = theAng(generator1); 
  double phi = phiAng(generator1);
  double xtran = distTran(generator1); 
  double ytran = distTran(generator1); 
  double ztran = distTran(generator1);
  
  point tranVec(xtran, ytran, ztran);
  point kv(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
  
  molCopyR.changeMoleculeMultiRotate(angle, kv, l, tranVec); // Adjust 'l' if needed
  bool cacaDist = molCopyR.checkCalphas(l); // Adjust 'l' if needed
  
  // if ca-ca distance restraint broken - return early
  if (cacaDist) { 
    
    return std::make_tuple(false, molCopyR, molFit, improvementIndex); 
  
  // if ca-ca distance restraint okay - move on
  } else {

      moleculeFitAndState molFitTmp = molFit;
      scatterFit = molFitTmp.getOverallFit(ed, mixtureList, helRatList, molCopyR, kmin, kmaxCurr, l); // Adjust 'l' if needed
      double uProb = distributionR(generator1);
      
      // if fit doesn't improve - return early
      if (!checkTransition(fitTemp.first, scatterFit.first, uProb, k, noScatterFitSteps)) {

        return std::make_tuple(false, molCopyR, molFit, improvementIndex); 

      // if fit improves - write files + log + update true
      } else {


          // Update mol, molFit, and other related states as needed
          // Log or output details as required

      }
  }



// referenced version! 

  // std::tuple<bool, ktlMolecule, std::pair<double, double>, moleculeFitAndState>
bool
performAffineTransformation(std::vector<ktlMolecule>& mol,
                            int sectionIndex,
                            const std::default_random_engine& generator1, 
                            const std::uniform_real_distribution<>& distTran, 
                            const std::uniform_real_distribution<>& rotAng, 
                            const std::uniform_real_distribution<>& theAng, 
                            const std::uniform_real_distribution<>& phiAng, 
                            const std::uniform_real_distribution<double>& distributionR, 
                            const std::default_random_engine& generator, 
                            moleculeFitAndState& molFit, 
                            const std::vector<double>& ed, 
                            const std::vector<double>& mixtureList, 
                            const std::vector<double>& helRatList, 
                            double kmin, 
                            double kmaxCurr,
                            int l,
                            int k,
                            double& scatterFit,
                            int noScatterFitSteps,
                            int& improvementIndex) {

  bool updated = false
  
  ktlMolecule molCopyR = mol[l];

  // all the ingredients for making a random rigid body transformation to our molecule!
  double angle = rotAng(generator1); 
  double theta = theAng(generator1); 
  double phi = phiAng(generator1);

  double xtran = distTran(generator1); 
  double ytran = distTran(generator1); 
  double ztran = distTran(generator1);
  
  point tranVec(xtran, ytran, ztran);
  point kv(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
  
  molCopyR.changeMoleculeMultiRotate(angle, kv, l, tranVec);
  bool cacaDist = molCopyR.checkCalphas(l);

  // if ca-ca distance restraint okay - move on!
  if (!cacaDist) { 

    moleculeFitAndState molFitTmp = molFit;
    std::pair<double,double> fitTemp = molFitTmp.getOverallFit(ed, mixtureList, helRatList, molCopyR, kmin, kmaxCurr, l); // Adjust 'l' if needed
    double uProb = distributionR(generator1);
      
    // if fit improves - write files + log + update true
    if (checkTransition(fitTemp.first, scatterFit.first, uProb, k, noScatterFitSteps)) {
        
      updated = true;
      improvementIndex++;

      scatterFit = fitTemp;
      mol[l] = molCopyR;
      molFit = molFitTmp;

      // Update mol, molFit, and other related states as needed
      // Log or output details as required

      }

  return updated;
  }

  

  // Instead of returning and reassigning - reference arguments
  // std::vector<ktlMolecule>& mol,  // Pass by reference to allow modifications
  }



// readin
  
    // Check if restarting from a previous run
    // if(strcmp(argv[3], "True") == 0) {

    //     std::stringstream filestring(argv[15]);
    //     std::string segment;
    //     std::vector<std::string> seglist;

    //     while(std::getline(filestring, segment, '+')) {
    //         seglist.push_back(segment);
    //     }

    //     for(size_t i = 0; i < seglist.size(); i++) {

    //         ktlMolecule molTmp;
    //         std::string sequenceLoc = std::string(argv[2]) + "fingerPrint" + std::to_string(i + 1) + ".dat";
    //         molTmp.readInSequence(sequenceLoc.c_str(), rmin, rmax, lmin);
    //         molTmp.readInCoordinates(seglist[i].c_str());
    //         molTmp.getHydrophobicResidues();
    //         mol.push_back(molTmp);
    //     }

    // // Fresh start
    // } else {

    //     int noStructures = std::atoi(argv[6]);
      
    //     for(int i = 0; i < noStructures; i++) {

    //         ktlMolecule molTmp;
    //         std::string sequenceLoc = std::string(argv[2]) + "fingerPrint" + std::to_string(i + 1) + ".dat";
    //         molTmp.readInSequence(sequenceLoc.c_str(), rmin, rmax, lmin);
    //         std::string coordinateLoc = std::string(argv[2]) + "coordinates" + std::to_string(i + 1) + ".dat";
    //         molTmp.readInCoordinates(coordinateLoc.c_str());
    //         molTmp.getHydrophobicResidues();
    //         mol.push_back(molTmp);
    //     }
    // }