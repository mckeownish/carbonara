#ifndef HYD_SH_M
#define HYD_SH_M

#include "ktlMoleculeRandom.h"
#include <algorithm>  
#include <functional> 


class hydrationShellMinimal{
public:
  /**
     Empty constructor
   **/   
  hydrationShellMinimal(){};
  /**
     @overload
     Constructor
     @param molIn existing molecule for which we create the hydration shell
     @param RinInp inner radius of hydration "shell"
     @param RoutInp outer radius of hydration "shell"
     @param RhyIn hydaration shell radius below which no sovent can appear
     @param ntrivsIn, the number of solvents in a ring surroudning the section of protein (6 as standard)
     @param helixRatioIn the number of rings of solvents as a fraction of the number of amino acids in the helical section
     @param solventsPerLinkIn, the number of rings per linker section (basically always 1)
     @param mutualDistCutOffIn the distance between which solvents where we class then as being in a shared space, we then choose only one of them.
     @param rmin min calpha distance
     @param rmax max calpha disatance
     @param lmin is the distacne where by two moelcules are considred too close
  **/  
  hydrationShellMinimal(ktlMolecule &molIn,double RinInp,double RoutInp,double &RhyIn,int ntrivsIn,double helixRatioIn,int solventsPerLinkIn,double &mutualDistCutOffIn,double &rmin,double &rmax,double &lmin);
  
  void getPointAndMidlengthMulti(int i,int &hIndex);
  void getPointAndMidlengthStraight(int &sec,int &part,int &hIndex,std::string soe,int &lenSec);
  void tubeParamList();
  void getAllHelices();
  point getDirec(int i);
  point getHydTan(int i);
  point getHydNorm(int i);
  point getHydBinorm(int i);
  double hydTubeLength(int i);
  point hydTubeCentre(int i);
  point getTangent(int index,int subindex);
  point getNormal(int index,int subindex);
  point getBinormal(int index,int subindex);
  point getCentrePoint(int index,int subIndex);
  void resetMolecule();
  int getMolSize();
  int getNoSections();
  int getNoKvals();
  void makeInitialSegData(point &cp,point &T,point &N1,double &tm,int index,int &nseg);
  void constructInitialState();
  void generateHydrationLayer();
  std::vector<std::vector<point> > getHydrationLayer(int i);
  void solventMoleculeDistances(std::vector<double> &molSolDistances,std::vector<double> &solSolDistances);
   std::vector<std::vector<point> > returnFlatSolList();
   void writeHydrationShellToFile(const char* filename);
 private:
  ktlMolecule mol;
  std::vector<point> direcList;
  std::vector<point> frameTan;
  std::vector<point> frameNorm;
  std::vector<point> frameBinorm;
  std::vector<point> midPointList;
  std::vector<std::string> nameLst;
  std::vector<int> lengthSec;
  std::vector<double> halfLengthList,cstepList;
  int molsize;
  double Pi,Rin,Rout,Rhy;
  std::vector<std::vector<std::vector<point> > > allSegments;
  std::vector<std::vector<std::vector<int> > > allTruthTables;
  int solventsPerLink;
  double helixRatio;
  int ntrivs;
  std::vector<std::vector<point> > helixpts;
  std::vector<std::vector<point> > helPtsFlat;
  std::vector<std::vector<point> > solPtsFlat;
  std::vector<std::vector<double> > distanceSets;
  std::vector<std::pair<double,double> > distanceList;
  int solIntbound,molIntbound,solmolIntbound,nscat,nSize;
  double maxDistChange,kminVal,kmaxVal,maxDistVal,maxdist,maxScatDist,mutualDistCutOff;
  bool isTooClose;
};

#endif
