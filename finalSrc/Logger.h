#ifndef LOGGER
#define LOGGER

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "parameters.h"
#include "moleculeFitAndState.h"


class Logger {

private:

    std::ofstream logFile;
    std::string escapeJSON(const std::string& s);

public:

    Logger(const std::string& filePath);
    ~Logger();

    void logEntry(const int& improvementIndex, const int& fitStep, const double& scatterFitFirst, 
                  const double& writhePenalty, const double& overlapPenalty, const double& distanceConstraints, 
                  const long& durationCount, const double& kmaxCurr, const std::string& scatterPath,
                  const std::string& moleculePath);

    void logMetadata(const std::string& run, ModelParameters params, const std::vector<double>& helRatList);

    void consoleInitial(const double& scatterFitFirst, const double& writhePenalty,
                                const double& overlapPenalty, const double& distanceConstraints);

    void consoleCurrentStep(int step, int index, double currFit);

    void consoleFitAttempt(int step, int improveIndex, ModelParameters params, double scatterFitFirst, double scatterFitSecond);



};

#endif