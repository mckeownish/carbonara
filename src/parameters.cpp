#include "parameters.h"

ModelParameters loadParameters(const char* argv[]) {

    ModelParameters params;

    // command-line arguments for scatter range
    params.kmin = std::atof(argv[9]);
    params.kmax = std::atof(argv[10]);

//    if (params.kmin < params.q_lim && params.kmax > params.q_lim) {
//        params.kmaxCurr = 0.15;
//    } else {
//        params.kmaxCurr = params.kmax;
//    }

    params.kmaxCurr = 0.25;

    // If resuming from a previous run, adjust parameters accordingly.
    if (strcmp(argv[3], "True") == 0) {
        params.kmaxCurr = std::atof(argv[24]);  // Ensure this is the correct argument index for kmaxCurr.
        std::cout << "kmax curr is " << params.kmaxCurr << "\n";
        // params.improvementIndex = std::atoi(argv[17]);  // Ensure this is the correct argument index for improvementIndex.
    }

    // command-line arguments for max number of steps
    params.noScatterFitSteps = std::atoi(argv[11]);

    // command-line arguments for rigid body transformations
    if(strcmp(argv[18],"True") == 0){
        params.affineTrans=true;
    } else {
        params.affineTrans=false;
    }

    // logic for helRatList --change?
    params.helRatList.push_back(0.5);

    params.basePath = argv[12];

    return params;
}
