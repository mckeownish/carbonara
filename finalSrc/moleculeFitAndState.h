#ifndef ALG_ROUTE
#define ALG_ROUTE

#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "writheFP.h"
#include "experimentalData.h"
#include "parameters.h"

class moleculeFitAndState{
public:
  moleculeFitAndState(std::vector<ktlMolecule> &mol, ModelParameters& params);
  std::vector<ktlMolecule> getMolecule();
  void updateMolecule(std::vector<ktlMolecule> &molNew);
  double calculateScattering(experimentalData &ed,double &kmin,double &kmax,std::vector<double> &mixtureVals);
  void writeScatteringToFile(experimentalData &ed,double &kmin,double &kmax,const char* filename);
  void writeHyrdationShellToFile(const char* filename,int &i);
  double getOverlapPenalty(double &closestApproachDist,std::vector<double> &overlapDists);
  double applyOverlapPenalty();
  double applyDistanceConstraints();
  double applyDistanceConstraints(ktlMolecule &molNew,int &i);
  void calculateMoleculeDistances(ktlMolecule &molNew,int &i);
  void calcuateHydrationDistances(hydrationShellMinimal &hs,int &i);
  void applyWritheConstraint();
  void calculateConnectionPenalty(ktlMolecule &molNew, int &chInd);
  std::pair<double,double>  getFit();
  double getWrithePenalty();
  double getOverlapPenalty();
  double getDistanceConstraints();
  void alterWritheSet(ktlMolecule &molNew,int &i);
  std::pair<double,double> getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,double &kmin,double &kmax);
  std::pair<double,double>  getOverallFitForceConnection(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,double &kmin,double &kmax);
  std::pair<double,double> getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
  std::pair<double,double>  getOverallFitForceConnection(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
  double currFit;
private:
  std::vector<std::vector<double> > molDists;
  std::vector<std::vector<double> > solDists;
  std::vector<std::vector<double> > solMolDists;
  std::vector<std::vector<double> > overlapDistSet;
  std::vector<int> molSize;
  std::vector<int> noSol;
  double maxDist;
  double hydroPhobicPacking;
  std::vector<std::vector<double> > originalWrithes;
  std::vector<std::vector<double> > currWrithes;
  std::vector<double> maxDistMol;
  std::vector<double> maxDistSol;
  std::vector<double> contactPredPen;
  double writhePenalty;
  double originalOverlapPenalty;
  double scatterAndHydrationConstraint;
  double Rin,Rout,RShell,ntrivs,closestApproachDist;
  double solventsPerLink,rmin,rmax,lmin;
  std::vector<double> percentageCombinations;
  std::vector<ktlMolecule> mol;
  double connectionPenalty;
  std::vector<double> connectionPenaltySet;
  std::vector<hydrationShellMinimal> hydrationShellBest;
};

#endif
