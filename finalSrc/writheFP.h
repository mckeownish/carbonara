#ifndef WRITHEFP_H
#define WRITHEFP_H

#include "point.h"
#include "ktlMoleculeRandom.h"
#include "skmt.h"

class writheFP{
public:
  writheFP();
  point dbold(std::vector<point>& pointList,int size,int i,int j);
  point dboldLk(std::vector<point>& pointList1,std::vector<point>& pointList2,int size1,int size2,int i,int j);
  point te(std::vector<point>& pointList,int size,int i);
  double mu(std::vector<point>& pointList,int size,int i,int k,int m, int j);
  double wij(std::vector<point>& pointList,int size,int i,int j);
  double lkij(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i,int j);
  double acosC(double temp);
  double DI(std::vector<point>& pointList);
  double DIAbs(std::vector<point>& pointList);
  //std::vector<double>  DI(std::vector<point>& pointList);
  double DIClosed(std::vector<point>& pointList);
  double muLink(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i, int j);
  double DIClosedLk(std::vector<point>& pointList,std::vector<point>& pointList2);
  std::vector<std::pair<std::pair<int,int>,double> > DIDownSample(std::vector<std::vector<point> >& pointListIn);
  std::vector<std::pair<std::pair<int,int>,double> > DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn);
  double DIDownSampleSingle(std::vector<std::vector<point> >& pointListIn);
  double DIDownSampleAbsSingle(std::vector<std::vector<point> >& pointListIn);
  double DIDownSampleAbsSingleSKMT(ktlMolecule&molin);
  void  DIDownSampleWrite(std::vector<std::vector<point> >& pointListIn,const char* filename);
  void  DIDownSampleAbsWrite(std::vector<std::vector<point> >& pointListIn,const char* filename);
private:
  double den,lambda;
  double writhepl;
  std::vector<point> dlist;
  std::vector<point> tanlist;
  std::vector<point> tderiv;
  std::vector<double> densityList;
  std::vector<int> turnlist;
  std::vector<int> turnPts;
  std::vector<double> arcLengthLst;
};

#endif
