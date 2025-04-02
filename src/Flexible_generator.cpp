#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h>
#include <vector>
#include "writheFP.h"


double calculateWrithePenalty(ktlMolecule &mol) {
    double writhePenalty = 0.0;
    writheFP wfp;
    
    // Loop through each chain in the molecule
    for(int chainNum = 1; chainNum <= mol.noChains(); chainNum++) {
        // Get coordinates for this chain
        std::vector<std::vector<point> > coords = mol.getSubsecCoordinates(chainNum);
        
        // Calculate writhe for this chain
        double chainWrithe = wfp.DIDownSampleAbsSingle(coords);
        
        // Get the length of this subsection
        double secLen = double(mol.getSubsecSize(chainNum));
        
        // Calculate lower bound - same formula as in original implementation
        double lowerBound = std::pow((secLen/7.5), 1.6) - 3.0;
        
        // Calculate penalty using sigmoid function
        // Same as original: 1.0/(1.0 + exp(20.0*(newWrithe-lowerBound)))
        writhePenalty += 1.0/(1.0 + std::exp(20.0*(chainWrithe - lowerBound)));
    }
    
    return writhePenalty;
}


/***********************************************
  argv[1] sequence file
  argv[2] initial prediction coords  file g
  argv[3] prediction output directory
  argv[4] change single section...
 **********************************************/


int main( int argc, const char* argv[] )
{  
    /*************************************

    set up model parameters 

    *************************************/

    double lmin=4.0; // closest distance two non adjactent local (same secondary unit) moelcules can get
    double rmin=3.7;double rmax=3.9; // max and min calpha-Dists
    double closestApproachDist=3.9; // closest distance two non adjactent non local moelcules (different secondary unit) can get

    /*************************************
     
    Generate new structures - changing only one section 

    *************************************/


    // initialise the molecule class
    ktlMolecule mol;

    // read in from sequence
    mol.readInSequence(argv[1],rmin,rmax,lmin);
        
    // read in coordinates
    mol.readInCoordinates(argv[2]);

    // Calculate initial writhe penalty
    double initialWrithePenalty = calculateWrithePenalty(mol);
    std::cout << "Initial writhe: " << initialWrithePenalty;

   for (int repeat=0; repeat<25; repeat++){
        
        // Make a copy of our precious
        ktlMolecule molcp = mol;
 
        // For now change section user input
        int change_section = std::atoi(argv[4]);
        //check which chain we are in
        int totNoSecs=mol.getSubsecSize(1);
        int chain_num=1;
        while(change_section>(totNoSecs-1)){
          chain_num++;
          //std::cout<<mol.getSubsecSize(chain_num)<<"\n";
          totNoSecs = totNoSecs+mol.getSubsecSize(chain_num);
        }
        totNoSecs=totNoSecs-mol.getSubsecSize(chain_num)-1;
        change_section = change_section-totNoSecs-1;
        molcp.changeMoleculeSingleMulti(change_section,chain_num);
 
        // c alphas adjacent dist check! 
        int sect_num = 1;
        bool cacaDist= molcp.checkCalphas(chain_num);

        double newWrithePenalty = calculateWrithePenalty(molcp);
        std::cout << "Updated writhe: " << newWrithePenalty;

        //std::cout<<cacaDist<<" "<<std::atoi(argv[4])<<" "<<chain_num<<"\n";
        if(cacaDist==false){

        // write the new coordinates!
        int repeatnum=repeat+1;
        
        // add the repeat number to stringstream
        std::stringstream ss;
        ss<<repeat;

        // make an output name for each repeat
        char outputname[100];
        strcpy(outputname,argv[3]);
        strcat(outputname,"_");
        std::string tempStr = ss.str();
        const char* str = tempStr.c_str(); 
        strcat(outputname,str);
        strcat(outputname,".dat");
        molcp.writeMoleculeToFile(outputname);
        
        }
   }
}
