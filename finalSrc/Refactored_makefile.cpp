#include "experimentalData.h"
#include "hydrationShellRandom.h"
#include "ktlMoleculeRandom.h"
#include "moleculeFitAndState.h"
#include <chrono>
#include <cstring>
#include <string.h>

using namespace std::chrono;

/***********************************************
Command-line arguments:

  argv[1] scattering data file
  argv[2] sequence file location
  argv[3] restart tag (use to start from existing prediction)
  argv[4] paired distances file (can be empty)
  argv[5] fixed sections file (again can be empty)
  argv[6] number of structures
  argv[7] request to apply hydrophobic covering WITHIN monomers will be a list
of sections on which to apply it. Will say none if not. -- Currently not used
  argv[8] request to apply hydrophobic covering BETWEEN monomers will be a list
of pairs to try to hydropobically pair. Will say none if not. -- currently not
used argv[9] kmin argv[10] kmax argv[11] Max number of fitting steps argv[12]
prediction file argv[13] scattering output file argv[14] mixture list file, a
list of sets of numbers indicatig the allowed set of mixture percentages of each
species (e.g. dimer 20 monomer 80) argv[15] previous fit string in form
fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat argv[16] log file
location argv[17] last line of the previous fit log, this is only used for a
restart if argv[3] = True argv[18] True if we want to apply affine rotations,
false otherwise.
**********************************************/

// Function to check if a transition between states should be accepted based on
// the simulated annealing criterion
bool checkTransition(double &chiSqVal, double &chiSqCurr, double &uniformProb,
                     int index, int &maxSteps) {
  // Simplified decision logic based on chi-squared values
  if (chiSqVal < chiSqCurr) {
    return true;
  } else {
    return false;
  }
}

// Function to sort a vector of moleculeFitAndState objects based on their
// current fit
void sortVec(std::vector<moleculeFitAndState> &mfs) {
  std::sort(mfs.begin(), mfs.end(),
            [](const moleculeFitAndState &x, const moleculeFitAndState &y) {
              return x.currFit < y.currFit;
            });
}

// Function to tokenize a string based on a delimiter, used for parsing input
// arguments
void tokenize(std::string &str, const char delim,
              std::vector<std::string> &out) {
  // construct a stream from the string
  std::stringstream ss(str);

  std::string s;
  while (std::getline(ss, s, delim)) {
    out.push_back(s);
  }
}

// Function to calculate the hydrophobic packing penalty
double getHydrophobicPackingPenalty(double &packValue) {
  return 0.00001 * std::exp(3.0 * (packValue - 1.6));
}

// Return list of ktlMolecules
std::vector<ktlMolecule> initializeMolecules(const char *argv[], double rmin,
                                             double rmax, double lmin) {

  std::vector<ktlMolecule> molecules;
  int noStructures = std::atoi(argv[6]);
  std::string sequenceFileLocation = argv[2];

  if (strcmp(argv[3], "True") == 0) {
    // Restart from a previous fit: parse the previous fit string and initialize
    // molecules accordingly
    std::vector<std::string> seglist;
    tokenize(argv[15], '+', seglist);
    for (int i = 0; i < seglist.size(); i++) {
      ktlMolecule molTmp;
      std::string sequenceLoc =
          sequenceFileLocation + "fingerPrint" + std::to_string(i + 1) + ".dat";
      molTmp.readInSequence(sequenceLoc.c_str(), rmin, rmax, lmin);
      molTmp.readInCoordinates(seglist[i].c_str());
      molTmp.getHydrophobicResidues();
      molecules.push_back(molTmp);
    }
  } else {
    // Regular start: initialize molecules based on sequence and structure files
    for (int i = 0; i < noStructures; i++) {
      ktlMolecule molTmp;
      std::string sequenceLoc =
          sequenceFileLocation + "fingerPrint" + std::to_string(i + 1) + ".dat";
      molTmp.readInSequence(sequenceLoc.c_str(), rmin, rmax, lmin);
      std::string coordinateLoc =
          sequenceFileLocation + "coordinates" + std::to_string(i + 1) + ".dat";
      molTmp.readInCoordinates(coordinateLoc.c_str());
      molTmp.getHydrophobicResidues();
      molecules.push_back(molTmp);
    }
  }
  return molecules;
}

// Returns list of fixed sections for each molecule
std::vector<std::vector<int>> readFixedSections(const std::string &baseLocation,
                                                int noStructures) {
  std::vector<std::vector<int>> fixedSecLists;

  for (int i = 0; i < noStructures; ++i) {
    std::vector<int> fixedSecList;
    std::string fixedSecLoc = baseLocation + "varyingSectionSecondary" +
                              std::to_string(i + 1) + ".dat";

    std::ifstream fixedSecFile(fixedSecLoc);
    std::string line;
    while (std::getline(fixedSecFile, line)) {
      std::stringstream ss(line);
      int index;
      while (ss >> index) {
        fixedSecList.push_back(index);
      }
    }

    if (!fixedSecList.empty()) {
      fixedSecLists.push_back(fixedSecList);
    } else {
      std::cerr << "Failed to open or read from file: " << fixedSecLoc
                << std::endl;
    }

    fixedSecFile.close();
  }

  return fixedSecLists;
}

// Function to load fixed distance constraints into molecules (if included)
void loadFixedDistanceConstraints(std::vector<ktlMolecule> &molecules,
                                  const std::string &baseLocation,
                                  bool applyConstraints) {
  if (!applyConstraints) {
    return; // Skip if constraints are not to be applied
  }

  for (size_t i = 0; i < molecules.size(); ++i) {
    std::string fixedDistLoc = baseLocation + "fixedDistanceConstraints" +
                               std::to_string(i + 1) + ".dat";
    // Attempt to load contact predictions for the molecule
    // Assuming ktlMolecule class has a method `loadContactPredictions` that
    // takes the file path as input
    molecules[i].loadContactPredictions(fixedDistLoc.c_str());
  }
}

// Function to read the permissible mixture list from a file
std::vector<std::vector<double>>
readPermissibleMixtureList(const std::string &filePath) {
  std::vector<std::vector<double>> mixtureList;
  std::ifstream permissibleMixtureFiles(filePath);
  std::string line;

  if (permissibleMixtureFiles.is_open()) {
    while (std::getline(permissibleMixtureFiles, line)) {
      std::vector<double> mixtureSet;
      std::stringstream ss(line);
      double index;
      while (ss >> index) {
        mixtureSet.push_back(index);
        // Skip any remaining whitespace in the stream to ensure correct parsing
        // of the next number.
        ss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
      }
      // Only add non-empty mixture sets to the list
      if (!mixtureSet.empty()) {
        mixtureList.push_back(mixtureSet);
      }
    }
    permissibleMixtureFiles.close();
  } else {
    std::cerr << "Failed to open mixture file: " << filePath << std::endl;
  }

  return mixtureList;
}

//
std::pair<std::pair<double, double>, double> calculateInitialScatterFit(
    double kmin, double kmax, moleculeFitAndState &molFit,
    const experimentalData &ed,
    const std::vector<std::vector<double>> &mixtureList,
    const std::vector<double> &helRatList,
    const std::string
        &predictionFileBasePath, // Base path for saving initial scatter context
    const char *restartTag,      // Indicates if this is a restart
    const char *restartKmax = nullptr,    // New kmax value if restarting
    const char *restartIndex = nullptr) { // Improvement index if restarting

  std::pair<double, double> scatterFit;
  double kmaxCurr;
  int improvementIndex = 0; // Initialize improvement index

  // Determine the appropriate kmaxCurr based on the condition
  if (kmin < 0.15 && kmax > 0.15) {
    kmaxCurr = 0.15;
  } else {
    kmaxCurr = kmax;
  }

  // Calculate the scatter fit using the determined kmaxCurr
  scatterFit =
      molFit.getOverallFit(ed, mixtureList, helRatList, kmin, kmaxCurr);

  // Construct the path for saving the initial scatter context
  std::string origMolLoc = predictionFileBasePath + "origScatter.dat";

  // Handle restart logic
  if (strcmp(restartTag, "True") == 0 && restartKmax != nullptr &&
      restartIndex != nullptr) {
    kmaxCurr = std::atof(restartKmax);
    improvementIndex = std::atoi(restartIndex);
    std::cout << "Restarting: kmax curr is " << kmaxCurr
              << ", improvement index is " << improvementIndex << "\n";
  }

  // Return both the scatter fit and the kmaxCurr used, along with the
  // improvement index
  return {{scatterFit, kmaxCurr}, improvementIndex};
}

/* -------------------------------------------------------------------------------------------
 */

int main(int argc, const char *argv[]) {

  // unsure
  bool doAll = false;

  // Base location for files
  std::string baseLocation = argv[2];

  // Model parameters setup
  double lmin = 4.0; // closest distance two non adjactent local (same secondary
                     // unit) moelcules can get
  double rmin = 3.7, rmax = 3.9; // Minimum and maximum C-alpha distances
  double closestApproachDist =
      3.9; // closest distance two non adjactent non local moelcules (different
           // secondary unit) can get

  // REFACTORED BITCHES

  // Initialize molecules
  std::vector<ktlMolecule> molecules =
      initializeMolecules(argv, rmin, rmax, lmin);

  int noStructures = std::atoi(argv[6]);

  // Call the readFixedSections function
  std::vector<std::vector<int>> fixedSecLists =
      readFixedSections(baseLocation, noStructures);

  // Check if fixed distance constraints should be applied
  bool applyConstraints = strcmp(argv[4], "True") == 0;

  // Load fixed distance constraints (contact predictions / sulfide bonds), if
  // applicable
  loadFixedDistanceConstraints(molecules, baseLocation, applyConstraints);

  // Read in the scattering and set up the scattering model
  experimentalData ed(argv[1]);

  std::string mixtureListFilePath = argv[14];

  // Load the permissible mixture list
  std::vector<std::vector<double>> mixtureList =
      readPermissibleMixtureList(mixtureListFilePath);

  // Initialise hydration shell parameters
  double Rin = 6.0, Rout = 7.0, RShell = 5.5;
  int ntrivs = 6, int solventsPerLink = 1;
  std::vector<double> helRatList = {0.5};
  double kmin = std::atof(argv[9]), kmax = std::atof(argv[10]);

  // calculate the initial fit, this includes all constraints
  moleculeFitAndState molFit(mol, Rin, Rout, RShell, ntrivs,
                             closestApproachDist, solventsPerLink, rmin, rmax,
                             lmin);

  // begin with up to 0.1
  std::pair<double, double> scatterFit;
  double kmaxCurr;
  if (kmin < 0.15 && kmax > 0.15) {
    kmaxCurr = 0.15;
    scatterFit =
        molFit.getOverallFit(ed, mixtureList, helRatList, kmin, kmaxCurr);
  } else {
    kmaxCurr = kmax;
    scatterFit = molFit.getOverallFit(ed, mixtureList, helRatList, kmin, kmax);
  }

  // save initial scatter for context

  char origMolLoc[100];
  strcpy(origMolLoc, argv[12]);
  strcat(origMolLoc, "origScatter.dat");

  int improvementIndex = 0;

  if ((strcmp(argv[3], "True") == 0)) {
    // hre we restart, we set kmaxCurr to be  the value previously obtained
    kmaxCurr = std::atof(argv[24]);
    std::cout << "kmax curr is " << kmaxCurr << "\n";
    improvementIndex = std::atoi(argv[17]);
  }

  // regenrate molecule hydration layer to update tha fit
  moleculeFitAndState molFitOg(mol, Rin, Rout, RShell, ntrivs,
                               closestApproachDist, solventsPerLink, rmin, rmax,
                               lmin);
  std::pair<double, double> scatterFitOut =
      molFitOg.getOverallFit(ed, mixtureList, helRatList, kmin, kmaxCurr);
  molFitOg.writeScatteringToFile(ed, kmin, kmaxCurr, origMolLoc);

  /**********************************************************************

   initialise the log file

  **********************************************************************/

  std::ofstream logFile;

  logFile.open(argv[16], std::ios::app);
  auto start = high_resolution_clock::now();
  auto curr = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(curr - start);

  if ((strcmp(argv[3], "True") == 1)) {
    // write the fit parameters if this is the first run
    logFile << "Run: " << argv[16] << " \n";
    logFile << "total q range considered " << kmin << " " << kmax << "\n";
    logFile << "initial q range considered " << kmin << " " << 0.2 << "\n";
    for (int i = 0; i < helRatList.size(); i++) {
      logFile << "helix solvent ratio " << i << " " << helRatList[i] << " ";
    }
    logFile << "\n";
    logFile << "\n\n\n";
    // write initial prediction
    logFile << 0 << " " << 0 << " " << scatterFit.first << " "
            << molFit.getWrithePenalty() << " " << molFit.getOverlapPenalty()
            << " " << molFit.getDistanceConstraints() << " " << duration.count()
            << " " << kmaxCurr << " " << origMolLoc << "\n";
    logFile.flush();
  } else {
    // here we are restarting we set kmax to
  }
  /****************************************************************************

     Main algorithm

  ***************************************************************************/
  //

  // //set up loop parameters
  int k = 0;

  // /* the vector noSections tells us how many subsections are in each moelcule
  //    e.g. for a monomer/dimer mixture noSections[0]=1,noSections[1]=2.
  //  */

  std::vector<int> noSections;
  for (int i = 0; i < mol.size(); i++) {
    int noSectionsTmp = mol[i].noChains();
    noSections.push_back(noSectionsTmp);
  }

  int noScatterFitSteps = std::atoi(argv[11]);

  // set up the vector of existing structures

  // initialise to the initial guess

  int noHistoricalFits = 1;

  std::vector<moleculeFitAndState> molFitAndStateSet;

  for (int i = 0; i < noHistoricalFits; i++) {
    molFitAndStateSet.push_back(molFit);
  }

  std::random_device rdev{};
  std::default_random_engine generator{rdev()};
  std::random_device rdev2{};
  std::default_random_engine generator1{rdev2()};

  // these are for the rotation/translations
  std::uniform_real_distribution<> distTran(-4.0, 4.0);
  std::uniform_real_distribution<> rotAng(0.0, 2.0);
  std::uniform_real_distribution<> theAng(0.0, 3.14159265359);
  std::uniform_real_distribution<> phiAng(0.0, 6.28318530718);

  bool affineTrans = false;
  if (strcmp(argv[18], "True") == 0) {
    affineTrans = true;
  }

  int improvementIndexTest = 0;

  while (k < noScatterFitSteps) {
    std::cout << kmin << " " << kmaxCurr << " " << kmax << " "
              << scatterFit.first << " " << scatterFit.second << "\n";
    if (scatterFit.second < 0.0005 ||
        (improvementIndexTest > std::round(noScatterFitSteps / 5) &&
         scatterFit.second < 0.0007)) {
      // if we have achieved a sufficiently good fit include more data.
      kmaxCurr = kmaxCurr + 0.01;
      if (kmaxCurr > kmax) {
        kmaxCurr = kmax;
      }
      improvementIndexTest = 0;
      // generate a new first fit.
      scatterFit = molFitAndStateSet[0].getOverallFit(
          ed, mixtureList, helRatList, kmin, kmaxCurr);
    }
    improvementIndexTest++;
    // loop over the molecules (e.g monomer and dimer fit

    std::uniform_real_distribution<double> distributionR(0.0, 1.0);

    // choose which fit to
    double p = 0.7 - 0.6 * (k / noScatterFitSteps);
    std::binomial_distribution<> changeIndexProbability(noHistoricalFits - 1,
                                                        p);
    int index = changeIndexProbability(generator);
    molFit = molFitAndStateSet[index];
    mol = molFit.getMolecule();
    scatterFit = molFit.getFit();
    for (int l = 0; l < mol.size(); l++) {
      int netIndex = 0;

      // loop over the sections of the given molecule (i.e. if its a monomer
      // this loop is tivial, but not for a multimer

      for (int i = 1; i <= noSections[l]; i++) {

        // check if we rotate/translate
        if (affineTrans == true) {
          ktlMolecule molCopyR = mol[l];
          double angle = rotAng(generator1);
          double theta = theAng(generator1);
          double phi = phiAng(generator1);
          double xtran = distTran(generator1);
          double ytran = distTran(generator1);
          double ztran = distTran(generator1);
          point tranVec(xtran, ytran, ztran);
          point kv(std::sin(theta) * std::cos(phi),
                   std::sin(theta) * std::sin(phi), std::cos(theta));
          molCopyR.changeMoleculeMultiRotate(angle, kv, i, tranVec);
          bool cacaDist = molCopyR.checkCalphas(i);
          if (cacaDist == false) {

            // calculate the new fit for this
            moleculeFitAndState molFitTmp = molFit;
            // calculate all amino acid distances for changed molecule
            std::pair<double, double> fitTemp = molFitTmp.getOverallFit(
                ed, mixtureList, helRatList, molCopyR, kmin, kmaxCurr, l);
            // check if we have imporved
            // std::cout<<"how change ? "<<fitTemp<<" "<<scatterFit<<"\n";
            double uProb = distributionR(generator);
            if (checkTransition(fitTemp.first, scatterFit.first, uProb, k,
                                noScatterFitSteps)) {
              scatterFit = fitTemp;
              mol[l] = molCopyR;
              molFit = molFitTmp;
              // to output during fitting to "show the process"
              improvementIndex++;
              for (int ii = 0; ii < mol.size(); ii++) {
                int ind1 = ii + 1;
                std::string outputMolLoc(argv[12]);
                outputMolLoc += "Substep_";
                outputMolLoc += std::to_string(ind1);
                outputMolLoc += "_";
                outputMolLoc += std::to_string(improvementIndex);
                outputMolLoc += ".dat";
                mol[ii].writeMoleculeToFile(outputMolLoc.c_str());
              }
              std::string outputScatLoc(argv[12]);
              outputScatLoc += "SubstepScatter_";
              outputScatLoc += std::to_string(improvementIndex);
              outputScatLoc += ".dat";
              molFitTmp.writeScatteringToFile(ed, kmin, kmaxCurr,
                                              outputScatLoc.c_str());
              curr = high_resolution_clock::now();
              duration = duration_cast<microseconds>(curr - start);
              logFile << improvementIndex << " " << k << " " << scatterFit.first
                      << " " << molFitTmp.getWrithePenalty() << " "
                      << molFitTmp.getOverlapPenalty() << " "
                      << molFitTmp.getDistanceConstraints() << " "
                      << duration.count() << " " << kmaxCurr << " "
                      << outputScatLoc << "\n";
              logFile.flush();
            }
          }
        }

        // net index tells us how far we are through the whole moelcule
        if (i > 1) {
          netIndex = netIndex + mol[l].getSubsecSize(i - 1);
        }

        // Now loop over the secondary structures of the given unit or section

        for (int j = 0; j < mol[l].getSubsecSize(i) - 1; j++) {
          // std::cout<<" mol "<<l<<" sec "<<i<<" has this many sections
          // "<<mol[l].getSubsecSize(i)<<"\n";
          int totalIndex = netIndex + j;
          // in this if statement we check which secondary sections are being
          // changed
          if ((doAll == true) ||
              (std::find(fixedSecLists[l].begin(), fixedSecLists[l].end(),
                         totalIndex) != fixedSecLists[l].end())) {
            // print statement currently in to check what we are changing is
            // correct
            // std::cout<<" section "<<totalIndex<<" of unit "<<i<<" "<<" sub
            // set number "<<totalIndex-netIndex<<" being altered
            // "<<mol.getSubsecSize(i)<<"\n";

            // copy the molecule to change it and test if we do better

            // now we change the section
            ktlMolecule molCopy = mol[l];
            int indexCh = totalIndex - netIndex;
            molCopy.changeMoleculeSingleMulti(indexCh, i);
            // this (checkCalphas) checks if there haven't been any rouge
            // sections created (some occasional flaws in the procedure which
            // are to be ironed out
            bool cacaDist = molCopy.checkCalphas(i);
            if (cacaDist == false) {

              // calculate the new fit for this
              moleculeFitAndState molFitTmp = molFit;
              // calculate all amino acid distances for changed molecule
              std::pair<double, double> fitTemp = molFitTmp.getOverallFit(
                  ed, mixtureList, helRatList, molCopy, kmin, kmaxCurr, l);
              // check if we have imporved
              // std::cout<<"how change ? "<<fitTemp<<" "<<scatterFit<<"\n";
              double uProb = distributionR(generator);
              if (checkTransition(fitTemp.first, scatterFit.first, uProb, k,
                                  noScatterFitSteps)) {
                scatterFit = fitTemp;
                mol[l] = molCopy;
                molFit = molFitTmp;
                // to output during fitting to "show the process"
                improvementIndex++;
                for (int ii = 0; ii < mol.size(); ii++) {
                  std::stringstream ss;
                  std::stringstream ss1;
                  int ind1 = ii + 1;
                  ss << ind1;
                  ss1 << improvementIndex;
                  char outputMolLoc[100];
                  strcpy(outputMolLoc, argv[12]);
                  strcat(outputMolLoc, "Substep_");
                  const char *impStr = ss.str().c_str();
                  strcat(outputMolLoc, impStr);
                  strcat(outputMolLoc, "_");
                  const char *str = ss1.str().c_str();
                  strcat(outputMolLoc, str);
                  strcat(outputMolLoc, ".dat");
                  mol[ii].writeMoleculeToFile(outputMolLoc);
                }
                std::stringstream ss1;
                ss1 << improvementIndex;
                char outputMolLoc[100];
                strcpy(outputMolLoc, argv[12]);
                strcat(outputMolLoc, "SubstepScatter_");
                const char *impStr = ss1.str().c_str();
                strcat(outputMolLoc, impStr);
                strcat(outputMolLoc, ".dat");
                molFitTmp.writeScatteringToFile(ed, kmin, kmaxCurr,
                                                outputMolLoc);
                curr = high_resolution_clock::now();
                duration = duration_cast<microseconds>(curr - start);
                logFile << improvementIndex << " " << k << " "
                        << scatterFit.first << " "
                        << molFitTmp.getWrithePenalty() << " "
                        << molFitTmp.getOverlapPenalty() << " "
                        << molFitTmp.getDistanceConstraints() << " "
                        << duration.count() << " " << kmaxCurr << " "
                        << outputMolLoc << "\n";
                logFile.flush();
              }
            }
          }
        }
      }
    }
    molFitAndStateSet[index] = molFit;
    molFitAndStateSet[index].updateMolecule(mol);
    sortVec(molFitAndStateSet);
    for (int i = 0; i < noHistoricalFits; i++) {
      std::cout << "step " << k << " " << i << " "
                << molFitAndStateSet[i].currFit << "\n";
    }
    k++;
  }
  mol = molFitAndStateSet[0].getMolecule();
  for (int i = 0; i < mol.size(); i++) {
    std::stringstream ss;
    int ind1 = i + 1;
    ss << ind1;
    char outputMolLoc[100];
    strcpy(outputMolLoc, argv[12]);
    strcat(outputMolLoc, "_");
    const char *str = ss.str().c_str();
    strcat(outputMolLoc, str);
    strcat(outputMolLoc, ".dat");
    mol[i].writeMoleculeToFile(outputMolLoc);
  }
  // regenrate molecule hydration layer to update tha fit
  moleculeFitAndState molFitOut(mol, Rin, Rout, RShell, ntrivs,
                                closestApproachDist, solventsPerLink, rmin,
                                rmax, lmin);
  scatterFitOut =
      molFitOut.getOverallFit(ed, mixtureList, helRatList, kmin, kmaxCurr);
  molFitOut.writeScatteringToFile(ed, kmin, kmaxCurr, argv[13]);
  improvementIndex++;
  logFile << improvementIndex << " " << k << " " << scatterFit.first << " "
          << molFitOut.getWrithePenalty() << " "
          << molFitOut.getOverlapPenalty() << " "
          << molFitOut.getDistanceConstraints() << " " << duration.count()
          << " " << kmaxCurr << " " << argv[13] << "\n";
  logFile.flush();
  logFile.close();
}

/* The refactor */

// Initialising the molecule - both new starts + restarts
// Inputs
//      sequenceFileLocation
//      restartTag
//      prevFitString
//      noStructures

std::vector<ktlMolecule>
initializeMolecules(const std::string &sequenceFileLocation,
                    const std::string &restartTag,
                    const std::string &prevFitString, int noStructures) {

  std::vector<ktlMolecule> molecules;
  if (restartTag == "True") {
    std::vector<std::string> seglist;
    tokenize(prevFitString, '+', seglist);
    for (int i = 0; i < seglist.size(); ++i) {
      ktlMolecule molTmp;
      std::string sequenceLoc =
          sequenceFileLocation + "fingerPrint" + std::to_string(i + 1) + ".dat";
      molTmp.readInSequence(sequenceLoc.c_str(), rmin, rmax, lmin);
      molTmp.readInCoordinates(seglist[i].c_str());
      molTmp.getHydrophobicResidues();
      molecules.push_back(molTmp);
    }
  } else {
    for (int i = 0; i < noStructures; ++i) {
      ktlMolecule molTmp;
      std::string sequenceLoc =
          sequenceFileLocation + "fingerPrint" + std::to_string(i + 1) + ".dat";
      std::string coordinateLoc =
          sequenceFileLocation + "coordinates" + std::to_string(i + 1) + ".dat";
      molTmp.readInSequence(sequenceLoc.c_str(), rmin, rmax, lmin);
      molTmp.readInCoordinates(coordinateLoc.c_str());
      molTmp.getHydrophobicResidues();
      molecules.push_back(molTmp);
    }
  }
  return molecules;
}

// Read in fixed sections
std::vector<std::vector<int>> readFixedSections(const std::string &baseLocation,
                                                int noStructures) {
  std::vector<std::vector<int>> fixedSecLists;
  for (int i = 0; i < noStructures; ++i) {
    std::vector<int> fixedSecList;
    std::string fixedSecLoc = baseLocation + "varyingSectionSecondary" +
                              std::to_string(i + 1) + ".dat";
    std::ifstream fixedSecFile(fixedSecLoc);
    std::string line;
    while (std::getline(fixedSecFile, line)) {
      if (!line.empty()) {
        std::stringstream ss(line);
        int index;
        while (ss >> index) {
          fixedSecList.push_back(index);
        }
      }
    }
    fixedSecLists.push_back(fixedSecList);
  }
  return fixedSecLists;
}
