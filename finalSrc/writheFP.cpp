#include "writheFP.h"

const double PI  =3.141592653589793238463;

writheFP::writheFP(){};

point writheFP::dbold(std::vector<point>& pointList,int size,int i,int j){
  i = i%(size-1);
  j = j%(size-1);
  return pointList[i+1].dif(pointList[j+1]);
}

point writheFP::dboldLk(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i,int j){
  i = i%(size-1);
  j = j%(size2-1);
  return pointList[i+1].dif(pointList2[j+1]);
}

point writheFP::te(std::vector<point>& pointList,int size,int i){
  point temp = dbold(pointList,size,i,i-1);
  temp.normalise();
  return temp;
}



double writheFP::acosC(double temp){
  if(temp < -1.0){
    return PI;
  }else if(temp > 1.0){
    return 0;
  }else{
    return acos(temp);
  }
}

double writheFP::mu(std::vector<point>& pointList,int size,int i,int k,int m, int j){
  point ti= te(pointList,size,i);
  point tj=te(pointList,size,j);
  point dkm = dbold(pointList,size,k,m);
  point un = ti.cross(dkm);
  point deux = dkm.cross(tj);
  if( !un.isNonzero() || !deux.isNonzero()){
    return 0;
  }else{
    un.normalise();
    deux.normalise();
    double temp = un.dotprod(deux);
    point cp = ti.cross(tj);
    double signProd = dkm.dotprod(cp); 
    if(std::abs(signProd)<0.000000001){
      return 0.0;
    }else{
      if(signProd>= 0){
	return acosC(temp);
      }else{
	return (-1.0)*acosC(temp);
      }
    }
   }
}

double writheFP::muLink(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i, int j){
  int iplus = i+1;
  int iminus = i-1;
  if(i==0){
    iminus=size-1;
  }
  if(i==size-1){
    iplus=0;
  }
  int jplus = j+1;
  int jminus = j-1;
  if(j==0){
    jminus=size2-1;
  }
  if(j==size2-1){
    jplus=0;
  }
  point tiplus= pointList[iplus].dif(pointList[i]);
  point timinus= pointList[i].dif(pointList[iminus]);
  point tjplus= pointList2[jplus].dif(pointList2[j]);
  point tjminus= pointList2[j].dif(pointList2[jminus]);
  double ds1 =0.5*(tiplus.length()+timinus.length());
  double ds2 =0.5*(tjplus.length()+tjminus.length());
  timinus.normalise();
  tjminus.normalise();
  tiplus.normalise();
  tjplus.normalise();
  point ti = (timinus+tiplus)*0.5;
  point tj = (tjminus+tjplus)*0.5;
  point dkm =pointList[i].dif(pointList2[j]);
  double dkmlen = dkm.length();
  double dkmCubed =dkmlen*dkmlen*dkmlen;
  point un = ti.cross(tj);
  return(ds1*ds2*un.dotprod(dkm)/dkmCubed);
}




double writheFP::wij(std::vector<point>& pointList,int size,int i,int j){
  double temp =0;
  if(j>i){
    double t1 = mu(pointList,size,i,i-1,j-1,j);
    double t2 = mu(pointList,size,i,i,j-1,j);
    double t3 = mu(pointList,size,i,i-1,j,j);
    double t4 = mu(pointList,size,i,i,j,j);
    temp = t1-t2-t3+t4;
  }
  return temp;
}


double writheFP::DI(std::vector<point>& pointList){
  //pointList.push_back(pointList[0]);
  int listSize = pointList.size();
  double sigsum=0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +wij(pointList,listSize,i,j);
     }
  }
  return sigsum/(2*PI);
}

double writheFP::DIAbs(std::vector<point>& pointList){
  //pointList.push_back(pointList[0]);
  int listSize = pointList.size();
  double sigsum=0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +std::abs(wij(pointList,listSize,i,j));
     }
  }
  return sigsum/(2*PI);
}

double writheFP::DIClosed(std::vector<point>& pointList){
  pointList.push_back(pointList[0]);
  int listSize = pointList.size();
  double sigsum=0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +wij(pointList,listSize,i,j);
     }
  }
  return sigsum/(2*PI);
}





double writheFP::DIClosedLk(std::vector<point>& pointList,std::vector<point>& pointList2){
  int listSize = pointList.size();
  int listSize2 = pointList2.size();
  double sigsum=0.0;
  for(int i=0;i<listSize;i++){
    for(int j=0;j<listSize2;j++){
      sigsum = sigsum +muLink(pointList,pointList2,listSize,listSize2,i,j);
     }
  }
  return sigsum/(4*PI);
}


std::vector<std::pair<std::pair<int,int>,double> > writheFP::DIDownSample(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DI(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         std::pair<std::pair<int,int>,double> fppoint;
         fppoint.first=pr;fppoint.second=wrval;
         fingerPrintList.push_back(fppoint);
       }
    }   
  }
  return fingerPrintList;
}




std::vector<std::pair<std::pair<int,int>,double> > writheFP::DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DIAbs(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         std::pair<std::pair<int,int>,double> fppoint;
         fppoint.first=pr;fppoint.second=wrval;
         fingerPrintList.push_back(fppoint);
       }
    }   
  }
  return fingerPrintList;
}

double  writheFP::DIDownSampleSingle(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  double wrval = DI(mol);
  return wrval;
}

double writheFP::DIDownSampleAbsSingle(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  double wrval = DIAbs(mol);
  return wrval;
}

double writheFP::DIDownSampleAbsSingleSKMT(ktlMolecule&molin){
  // first downsample
  skmt s;
  std::vector<point> mol = s.getSKMTCurve(molin);
  double wrval = DIAbs(mol);
  return wrval;
}




void writheFP::DIDownSampleWrite(std::vector<std::vector<point> >& pointListIn,const char* filename){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  std::ofstream fingerFile;
  fingerFile.open(filename);
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DI(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         fingerFile<<k<<" "<<k+l<<" "<<wrval<<"\n";
       }
    }   
  }
  fingerFile.close();
}

void writheFP::DIDownSampleAbsWrite(std::vector<std::vector<point> >& pointListIn,const char* filename){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  std::ofstream fingerFile;
  fingerFile.open(filename);
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DIAbs(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         fingerFile<<k<<" "<<k+l<<" "<<wrval<<"\n";
       }
    }   
  }
  fingerFile.close();
}
