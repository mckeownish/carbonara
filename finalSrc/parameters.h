#ifndef PARAMS
#define PARAMS

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

struct ModelParameters {

    double lmin = 4.0;
    double rmin = 3.7;
    double rmax = 3.9;
    double closestApproachDist = 3.9;
    double Rin = 6.0;
    double Rout = 7.0;
    double RShell = 5.5;
    int ntrivs = 6;
    int solventsPerLink = 1;
    double q_lim = 0.15;

    // user definable
    double kmin = 0.0;
    double kmax = 0.0;
    double kmaxCurr = 0.0;
    int noScatterFitSteps = 5000;
    bool affineTrans=false;

    std::vector<double> helRatList;
    std::vector< std::vector<double> > mixtureList;
    std::string basePath;

    int improvementIndexTest = 0;
    int noHistoricalFits = 1;
    
};

// All model parameters
ModelParameters loadParameters(const char* argv[]);

#endif