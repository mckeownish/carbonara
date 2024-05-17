
#include "hydrationShellRandom.h"


/***********

Constructor, needs to be here, takes an existing moleclule (the ktlMolecule object which is passed by reference) the parameters RinInp,RoutInp,RhyIn,ntrivsIn,helixRatioIn,solventsPerLinkIn,mutualDistCutOffIn,rmin,rmax,lmin, are all specified in the main file and their values have been determined empirically in the main paper.

 *******/


hydrationShellMinimal::hydrationShellMinimal(ktlMolecule &molIn,double RinInp,double RoutInp,double &RhyIn,int ntrivsIn,double helixRatioIn,int solventsPerLinkIn,double &mutualDistCutOffIn,double &rmin,double &rmax,double &lmin){
  mol = molIn;
  molsize = mol.molSize();
  nSize = mol.noSecSize();
  direcList.resize(molsize);
  frameTan.resize(molsize);
  frameNorm.resize(molsize);
  frameBinorm.resize(molsize);
  halfLengthList.resize(molsize);
  nameLst.resize(molsize);
  lengthSec.resize(molsize);
  midPointList.resize(molsize);
  Pi = 3.1415926535897932384626;
  Rin = RinInp;Rout = RoutInp;Rhy =RhyIn;
  allSegments.resize(molsize);
  allTruthTables.resize(molsize);
  cstepList.resize(molsize);
  solventsPerLink = solventsPerLinkIn;
  helixRatio = helixRatioIn;
  ntrivs = ntrivsIn;
  helixpts.resize(molsize);
  //declare the empty distances list
  maxDistChange =0.0;
  mutualDistCutOff =mutualDistCutOffIn;
  isTooClose= false;
}


point hydrationShellMinimal::getCentrePoint(int index,int subIndex){
  return mol.getCoordinate(index,subIndex);
}



/**********************************

Takes a linker section (of index i) and finds its effective cylindrical shell on which to place the solvent molcules. h inedex is the index of the hdration shell (a single linker of index i has as many subsections as its number of amino acids in the hIndex).

 *********************************/

void hydrationShellMinimal::getPointAndMidlengthMulti(int i,int &hIndex){
  std::vector<point> coords  = mol.getCoordinatesSection(i);int sz =coords.size();
  std::vector<point> lines;point mean(0.0,0.0,0.0);
  for(int j=0;j<sz;j++){
    mean = mean + coords[j];
  }
  mean = mean*(1.0/double(sz));
  point endDif = mean - coords[sz-1];
  for(int j=0;j<sz;j++){
    lines.push_back(mean - coords[j]);
    if(lines[j].dotprod(endDif)<0.0){
      lines[j] = lines[j]*(-1.0);
    }
  }
  point direc(0.0,0.0,0.0);
  for(int j=0;j<sz;j++){
    direc = direc + lines[j];
  }
  double hlfDist= 0.5*(mean.eDist(coords[0])+mean.eDist(coords[sz-1]));
  direc.normalise();
  midPointList[hIndex]=mean;
  frameTan[hIndex] = direc;
  point norm;
  point vert(0.0,0.0,1.0);
  if(std::abs(direc.dotprod(vert))<0.99999999){
    norm = vert.cross(direc);
    norm.normalise();
  }else{
    point vert2(0.0,1.0,0.0);
    norm= vert2.cross(direc);
    norm.normalise();
  }
  frameNorm[hIndex] = norm;
  frameBinorm[hIndex] = direc.cross(norm);
  direcList[hIndex]=direc;
  halfLengthList[hIndex]=hlfDist;
}


void hydrationShellMinimal::getPointAndMidlengthStraight(int &sec,int &part,int &hIndex,std::string soe,int &lenSec){
  point startPt,endPt,midPoint;
  if(soe == "start"){
    if(sec==0){
      midPoint = mol.getCoordinate(sec,part);
      if(lenSec==1){
	endPt = mol.getCoordinate(sec+1,0);
      }else{
	endPt = mol.getCoordinate(sec,part+1);
      }
      startPt = midPoint -(endPt-midPoint);
    }else{
      if(lenSec==1){
	int noPrev =mol.getUnitNo(sec-1);
	startPt = mol.getCoordinate(sec-1,noPrev-1);
	midPoint = mol.getCoordinate(sec,part);
	endPt = midPoint+(midPoint-startPt);
      }else{
	int noPrev =mol.getUnitNo(sec-1);
	startPt = mol.getCoordinate(sec-1,noPrev-1);
	midPoint = mol.getCoordinate(sec,part);
	endPt = mol.getCoordinate(sec,part+1);
      }
    }
  }else if(soe == "end"){
    if(sec==nSize-1){
      if(lenSec==1){
	int noPrev =mol.getUnitNo(sec-1);
	startPt = mol.getCoordinate(sec-1,noPrev-1);
      }else{
	startPt = mol.getCoordinate(sec,-1);
      }
      midPoint = mol.getCoordinate(sec,part);
      endPt = midPoint+(midPoint-startPt);
    }else{
       if(lenSec==1){
	 startPt = midPoint -(endPt-midPoint);
	 midPoint = mol.getCoordinate(sec,part);
	 endPt = mol.getCoordinate(sec+1,0);
       }else{
	 startPt = mol.getCoordinate(sec,part-1);
	 midPoint = mol.getCoordinate(sec,part);
	 endPt = mol.getCoordinate(sec+1,0);
       }
    }
  }else{
      startPt = mol.getCoordinate(sec,part-1);
      midPoint = mol.getCoordinate(sec,part);
      endPt = mol.getCoordinate(sec,part+1);
  }
  double endDist = 0.5*startPt.eDist(endPt);
  point direc = endPt-startPt;
  direc.normalise();
  direcList[hIndex]=direc;
  frameTan[hIndex] = direc;
  double normt;
  if(sec==0||sec == nSize-1){
    point vert(0.0,0.0,1.0);
    if(std::abs(direc.dotprod(vert))<0.9999999){
      point N1 = vert.cross(direc);
      N1.normalise();
      frameNorm[hIndex] = N1;
    }else{
      point vert2(0.0,1.0,0.0);
      point N1 = vert2.cross(direc);
      N1.normalise();
      frameNorm[hIndex] = N1;
    }
  }else{
    point bmina = endPt-startPt;
    normt = (midPoint-startPt).dotprod(bmina);
    normt = normt/(bmina.dotprod(bmina));
    point frnm =  startPt + (endPt-startPt)*normt-midPoint;
    frnm.normalise();
    frameNorm[hIndex] =frnm;
  }
  frameNorm[hIndex].normalise();
  frameBinorm[hIndex] = frameTan[hIndex].cross(frameNorm[hIndex]);
  midPointList[hIndex]= startPt;   
  halfLengthList[hIndex]=endDist*0.5;
}


/**********************************

tubeParamList runs over the whole molecule and usese the shell creation algorithms to create the mathematical description of the backbone (basically descriptions of a bunch of cylinders surrounding the c-alpha backbone.

 *********************************/

void hydrationShellMinimal::tubeParamList(){
  int currHIndex=0;
  for(int i=0;i<nSize;i++){
    if(mol.getType(i)=="Helix" || mol.getType(i)=="Strand"){
      getPointAndMidlengthMulti(i,currHIndex);
      nameLst[currHIndex] =mol.getType(i);
      lengthSec[currHIndex] = mol.getUnitNo(i);
      currHIndex++;
    }else{
      // loop over the whole section
      for(int j=0;j<mol.getUnitNo(i);j++){
	nameLst[currHIndex] = "loop";
        lengthSec[currHIndex] = 1;
	if(j==0){
	  getPointAndMidlengthStraight(i,j,currHIndex,"start",lengthSec[currHIndex]);
	}else if(j==mol.getUnitNo(i)-1){
	  getPointAndMidlengthStraight(i,j,currHIndex,"end",lengthSec[currHIndex]);
	}else{
	  getPointAndMidlengthStraight(i,j,currHIndex,"middle",lengthSec[currHIndex]);
	}
	currHIndex++;
      }
    }
  }
}


/**********************************

tubeParamList runs over the whole molecule and usese the shell creation algorithms to create the full (initial) hydration shell

 *********************************/

void hydrationShellMinimal::getAllHelices(){
  int currHIndex=0;
  for(int i=0;i<nSize;i++){
    if(mol.getType(i)=="Helix" || mol.getType(i)=="Strand"){
      helixpts[currHIndex]=mol.getCoordinatesSection(i);
      currHIndex++;
    }else{
      // loop over the whole section
      for(int j=0;j<mol.getUnitNo(i);j++){
	std::vector<point> helpt;
      	helpt.push_back(mol.getCoordinate(i,j));
        helixpts[currHIndex] = helpt;
      	currHIndex++;
      }
    }
  }
};

// methods to call components of the hydration model

point hydrationShellMinimal::getDirec(int i){
  return direcList[i];
}

point hydrationShellMinimal::getHydTan(int i){
  return frameTan[i];
}

point hydrationShellMinimal::getHydNorm(int i){
  return frameNorm[i];
}

point hydrationShellMinimal::getHydBinorm(int i){
  return frameBinorm[i];
}

double hydrationShellMinimal::hydTubeLength(int i){
  return 2*halfLengthList[i];
}

point hydrationShellMinimal::hydTubeCentre(int i){
  return midPointList[i];
}

int hydrationShellMinimal::getMolSize(){
  return mol.molSize();
}

int hydrationShellMinimal::getNoSections(){
  return mol.noSecSize();
}


point hydrationShellMinimal::getTangent(int secindex,int subIndex){
  return mol.getTangent(secindex,subIndex);
};

point hydrationShellMinimal::getNormal(int secindex,int subIndex){
  return mol.getNormal(secindex,subIndex);
};

point hydrationShellMinimal::getBinormal(int secindex,int subIndex){
  return mol.getBinormal(secindex,subIndex);
};



/**********************************

tubeParamList runs over the whole molecule and usese the shell creation algorithms to create the full (initial) hydration shell

*********************************/

void hydrationShellMinimal::makeInitialSegData(point &cp,point &T,point &N1,double &tm,int index,int &nseg){
  point N12 = T.cross(N1);
  double cstep;
  if(nseg>1){
    cstep = 2*tm/double(nseg-1);
  }else{
    cstep =0.0;
  }
  //find the points along the "tube" axis where the hydrations discs will be centered
  cstepList[index] = cstep;
  std::vector<point> line;
  /*if(index<4){
    cp.printPoint();
    T.printPoint();
    N1.printPoint();
    }*/
  if(nseg%2==0){
    for(int i=-nseg/2;i <= -1;i++){
      line.push_back(cp + T*(i*cstep + 0.5*cstep));
    }
    for(int i=1;i<= nseg/2;i++){
      line.push_back(cp + T*(i*cstep - 0.5*cstep));
    }
  }else{
    if(nseg>1){
      for(int i=-(nseg-1)/2;i<=(nseg-1)/2;i++){
        line.push_back(cp + T*(i*cstep));
      }
    }else{
        line.push_back(cp);
    }
  }
  std::vector<std::vector<point> > segments;
  std::vector<std::vector<int> > truthTable(line.size(),std::vector<int>(ntrivs,1));
  point p;
  std::vector<int> checkList;
  checkList.push_back(1000000);
  std::vector<std::vector<int> > checkDisc;
  std::vector<std::vector<std::vector<int> > > checkFull;
  for(int i=0;i<line.size();i++){
    std::vector<point> seg(ntrivs);
    for(int j=1;j<=ntrivs;j++){
      p = line[i] - T*0.5*cstep + (N1*std::cos(2*Pi*j/ntrivs + Pi/ntrivs) + N12*std::sin(2*Pi*j/ntrivs + Pi/ntrivs))*Rhy;
      seg[j-1] =p;
      checkDisc.push_back(checkList);
    }
    checkFull.push_back(checkDisc);
    segments.push_back(seg);
  }
  allSegments[index]= segments;
  allTruthTables[index] = truthTable;
}

void hydrationShellMinimal::constructInitialState(){
  for(int i=0;i<molsize;i++){
    point outerNorm= frameNorm[i]*(-1.0);
    if(nameLst[i]=="Helix"||nameLst[i]=="Strand"){
      // helix section
      //int noSegs = std::max(1,int(ceil(helixRatio*noCalphas)));
      int noCalphas = lengthSec[i];
      int noSegs = std::max(1,int(round(helixRatio*noCalphas)));
      //std::cout<<i<<" "<<noSegs<<" "<<" "<<helixRatio<<" "<<noCalphas<<"\n";
      makeInitialSegData(midPointList[i],direcList[i],outerNorm,halfLengthList[i],i,noSegs);
    }else{
      // linker
      makeInitialSegData(midPointList[i],direcList[i],outerNorm,halfLengthList[i],i,solventsPerLink);
    }
  }
}



std::vector<std::vector<point> > hydrationShellMinimal::getHydrationLayer(int i){
  return allSegments[i];
}


void hydrationShellMinimal::generateHydrationLayer(){
  tubeParamList();
  constructInitialState();
  getAllHelices();
}




void hydrationShellMinimal::solventMoleculeDistances(std::vector<double> &molSolDistances,std::vector<double> &solSolDistances){
  // first find all solvent moelcule distances
  for(int k=0;k<allSegments.size();k++){
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	bool overlapped = false;
	std::vector<double> solMolPosDists;
	for(int l=0;l<helixpts.size();l++){
	  for(int m=0;m<helixpts[l].size();m++){
	    // declare molpt
	    /*if(k==19 && i ==0  && l ==238 && m==0){
	       allSegments[k][0][j].printPoint();
	       helixpts[l][m].printPoint();
	       std::cout<<"dist check  "<<helixpts[l][m].eDist(allSegments[k][i][j])<<"\n";
	       }*/
	    if(k!=l){
	      double dist =helixpts[l][m].eDist(allSegments[k][i][j]);
	      if(dist<5.4){
		overlapped=true;
		//noSegsTot++;
		//brake++;
		/*if(-1<brake && brake <10){  
		   std::cout<<dist<<" "<<k<<" "<<i<<" "<<j<<" "<<l<<" "<<m<<"\n";
		   }*/
	      }else{
		solMolPosDists.push_back(dist);
	      }
	    }else{
	       double dist =helixpts[l][m].eDist(allSegments[k][i][j]);
	       solMolPosDists.push_back(dist);
	    }
	  }
	}
	if(overlapped==true){
	  allTruthTables[k][i][j]=0;
	}else{
	  // this solvent is still in play, store its distances.
	  if(molSolDistances.size()==0){
	    molSolDistances= solMolPosDists;
	  }else{
	    molSolDistances.insert(molSolDistances.end(),solMolPosDists.begin(),solMolPosDists.end());
	  }
	}
	//std::cout<<i<<" "<<j<<" "<<k<<" "<<molSolDistances.size()<<"\n";
      }
    }
  }
  //std::cout<<" mol sol here ? "<<noSegsTot<<" "<<molSolDistances.size()<<"\n";
  for(int k=0;k<allSegments.size()-1;k++){
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	//std::cout<<allTruthTables[k][i][j]<<"\n";
	if(allTruthTables[k][i][j]==1){
	  for(int l=k+1;l<allSegments.size();l++){
	    for(int m=0;m<allSegments[l].size();m++){
	      for(int n=0;n<allSegments[l][0].size();n++){
		if(allTruthTables[l][m][n]==1){
		 double soldist = allSegments[l][m][n].eDist(allSegments[k][i][j]);
		 //std::cout<<soldist<<"\n";
		 if(soldist<2.5){
		   allTruthTables[l][m][n]=0;
		 }else{
		   solSolDistances.push_back(soldist);
		 }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //std::cout<<"escape ? "<<molSolDistances.size()<<" "<<solSolDistances.size()<<"\n";
  // now we check for solvent solvent distances *between* sections (not within a section)
}



std::vector<std::vector<point> > hydrationShellMinimal::returnFlatSolList(){
  solPtsFlat.clear();
  int noflps =0;
  for(int k=0;k<allSegments.size();k++){
    std::vector<point> spts;
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	//std::cout<<i<<" "<<j<<" "<<k<<" "<<allTruthTables[k][i][j]<<"\n";
	if(allTruthTables[k][i][j]==1){
	  spts.push_back(allSegments[k][i][j]);
	}
      }
    }
    solPtsFlat.push_back(spts);
    //noflps = noflps + spts.size();
    //std::cout<<"in full segment "<<k<<" "<<spts.size()<<"\n";
      //also add the helical molecules to the flat helix list
  }
  // std::cout<<"full calc "<<noflps<<"\n";
  return solPtsFlat;
}




void hydrationShellMinimal::writeHydrationShellToFile(const char* filename){
  std::ofstream ofl;
  ofl.open(filename);
  if(ofl.is_open()){
    for(int k=0;k<allSegments.size();k++){
      for(int i=0;i<allSegments[k].size();i++){
	for(int j=0;j<allSegments[k][0].size();j++){
	  if(allTruthTables[k][i][j]==1){
	    ofl<<allSegments[k][i][j].getX()<<" "<<allSegments[k][i][j].getY()<<" "<<allSegments[k][i][j].getZ()<<"\n";
	  }
	}
      }
    }
  }
}

