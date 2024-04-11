#include "Logger.h"

Logger::Logger(const std::string& filePath) {
    logFile.open(filePath, std::ios::out | std::ios::app);
    if (!logFile.is_open()) {
        std::cerr << "Failed to open log file: " << filePath << std::endl;
    }
}

Logger::~Logger() {
    if (logFile.is_open()) {
        logFile.close();
    }
}

void Logger::logEntry(const int& improvementIndex, const int& fitStep, const double& scatterFitFirst, 
                      const double& writhePenalty, const double& overlapPenalty, const double& distanceConstraints, 
                      const long& durationCount, const double& kmaxCurr, const std::string& scatterPath,
                      const std::string& moleculePath) {

    if (logFile.is_open()) {

        logFile << "{";
        logFile << "\"ImprovementIndex\": " << improvementIndex << ", ";
        logFile << "\"FitStep\": " << fitStep << ", ";
        logFile << "\"ScatterFitFirst\": " << scatterFitFirst << ", ";
        logFile << "\"WrithePenalty\": " << writhePenalty << ", ";
        logFile << "\"OverlapPenalty\": " << overlapPenalty << ", ";
        logFile << "\"DistanceConstraints\": " << distanceConstraints << ", ";
        logFile << "\"DurationCount\": " << durationCount << ", ";
        logFile << "\"KmaxCurr\": " << kmaxCurr << ", ";
        logFile << "\"ScatterPath\": \"" << escapeJSON(scatterPath) << "\"" << ", ";
        logFile << "\"MoleculePath\": \"" << escapeJSON(moleculePath) << "\"";
        logFile << "}\n";

        logFile.flush();
    }
}

void Logger::logMetadata(const std::string& run, ModelParameters params,
                         const std::vector<double>& helRatList) {

    if (logFile.is_open()) {

        logFile << "{";
        logFile << "\"Run\": \"" << escapeJSON(run) << "\", ";
        logFile << "\"qmin\": \"" << params.kmin << "\", ";
        logFile << "\"qmax\": \"" << params.kmax << "\", ";
        logFile << "\"qmax_curr\": \"" << params.kmaxCurr << "\", ";

        logFile << "\"HelixSolventRatios\": [";
        for (size_t i = 0; i < helRatList.size(); ++i) {
            logFile << helRatList[i];
            if (i < helRatList.size() - 1) {
                logFile << ", ";
            }
        }
        logFile << "], ";

        logFile << "\"lmin\":" << params.lmin << ", ";
        logFile << "\"rmin\":" << params.rmin << ", ";
        logFile << "\"rmax\":" << params.rmax << ", ";
        logFile << "\"closestApproachDist\":" << params.closestApproachDist << ", ";
        logFile << "\"lmin\":" << params.lmin << ", ";
        logFile << "\"Rin\":" << params.Rin << ", ";
        logFile << "\"Rout\":" << params.Rout << ", ";
        logFile << "\"RShell\":" << params.RShell << ", ";
        logFile << "\"ntrivs\":" << params.ntrivs << ", ";
        logFile << "\"solventsPerLink\":" << params.solventsPerLink << "}";
        logFile << "\n";

        logFile.flush();
    }
}

std::string Logger::escapeJSON(const std::string& s) {

    std::ostringstream o;

    for (auto c : s) {

        switch (c) {
            case '"': o << "\\\""; break;
            case '\\': o << "\\\\"; break;
            case '\b': o << "\\b"; break;
            case '\f': o << "\\f"; break;
            case '\n': o << "\\n"; break;
            case '\r': o << "\\r"; break;
            case '\t': o << "\\t"; break;
            default: o << c; break;
        }
    }

    return o.str();
}

void Logger::consoleInitial(const double& scatterFitFirst, const double& writhePenalty,
                            const double& overlapPenalty, const double& distanceConstraints) {
                                
    std::cout << std::left << std::setw(80) << "Initial molecule"
                           << "\n";

    std::cout << std::left << std::setw(20) << "Scattering Fit" 
                           << std::setw(20) << "Overlap Penalty"
                           << std::setw(20) << "Writhe Penalty"
                           << std::setw(20) << "Contact Penalty"
                           << "\n";

    std::cout << std::left << std::setw(20) << scatterFitFirst
                           << std::setw(20) << overlapPenalty
                           << std::setw(20) << writhePenalty
                           << std::setw(20) << distanceConstraints
                           << "\n \n";

}

void Logger::consoleCurrentStep(int step, int index, double currFit) {

    std::cout << std::left << std::setw(10) << "Step"
                            << std::setw(10) << step 
                            << std::setw(10) << "Index" 
                            << std::setw(10) << index 
                            << std::setw(15) << "Current Fit" 
                            << std::setw(10) << std::setprecision(5) << currFit 
                            << "\n";

}   

void Logger::consoleFitAttempt(int step, int improveIndex, ModelParameters params, double scatterFitFirst, double scatterFitSecond) {
    
    if (step==0){
        std::cout << std::left << std::setw(20) << "Improved Index" 
                               << std::setw(20) << "Current Fit Step"
                               << std::setw(20) << "Current Scatter"
                               << std::setw(20) << "Updated Scatter"
                               << "\n";

    }

    
    std::cout << std::left
                
                << std::setw(20)  << improveIndex
                << std::setw(20)  << step
                << std::setw(20) << scatterFitFirst
                << std::setw(20) << scatterFitSecond
                << "\n";
    }