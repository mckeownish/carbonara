#include "ktlMoleculeRandom.h"
#include <iostream>

#ifdef __APPLE__
    #define EXPORT __attribute__((visibility("default")))
#else
    #define EXPORT
#endif

extern "C" {
    // Structure to return coordinates
    struct CoordArray {
        double* data;
        int size;
    };

    static ktlMolecule* currentMolecule = nullptr;
    static ktlMolecule* tempMolecule = nullptr;  // For storing temporary changes

    EXPORT bool init_state(const char* seq_file, const char* coord_file) {
        try {
            if (currentMolecule != nullptr) {
                delete currentMolecule;
            }
            if (tempMolecule != nullptr) {
                delete tempMolecule;
            }
            
            currentMolecule = new ktlMolecule();
            
            double rmin = 3.7;
            double rmax = 3.9;
            double lmin = 4.0;
            
            currentMolecule->readInSequence(seq_file, rmin, rmax, lmin);
            currentMolecule->readInCoordinates(coord_file);
            
            // Initialize temp molecule
            tempMolecule = new ktlMolecule(*currentMolecule);
            
            return true;
        } catch (...) {
            return false;
        }
    }

//     EXPORT CoordArray try_update(int change_section) {
// //        std::cout << "Trying update..." << std::endl;
        
//         CoordArray result = {nullptr, 0};
//         if (!currentMolecule || !tempMolecule) return result;

//         // Reset temp molecule to current state
//         *tempMolecule = *currentMolecule;
        
//         // Calculate chain number
//         int totNoSecs = tempMolecule->getSubsecSize(1);
//         int chain_num = 1;
//         while(change_section > (totNoSecs-1)) {
//             chain_num++;
//             totNoSecs = totNoSecs + tempMolecule->getSubsecSize(chain_num);
//         }
//         totNoSecs = totNoSecs - tempMolecule->getSubsecSize(chain_num) - 1;
//         change_section = change_section - totNoSecs - 1;
        
//         // Try the change on temp molecule
//         tempMolecule->changeMoleculeSingleMulti(change_section, chain_num);
        
//         if (!tempMolecule->checkCalphas(chain_num)) {
//             // Get coordinates from successful change
//             auto coords = tempMolecule->getCoordinates();
            
//             // Calculate total size needed
//             int total_size = 0;
//             for (auto& chain : coords) {
//                 total_size += chain.size() * 3;
//             }
            
//             // Allocate array
//             double* data = new double[total_size];
//             int idx = 0;
            
//             // Copy coordinates to flat array
//             for (auto& chain : coords) {
//                 for (auto& point : chain) {
//                     data[idx++] = point.getX();
//                     data[idx++] = point.getY();
//                     data[idx++] = point.getZ();
//                 }
//             }
            
//             result.data = data;
//             result.size = total_size;
//         }
        
//         return result;
//     }

EXPORT CoordArray try_update(int change_section) {
    CoordArray result = {nullptr, 0};
    if (!currentMolecule || !tempMolecule) return result;

    const int MAX_ATTEMPTS = 100;
    for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
        // Reset temp molecule to current state
        *tempMolecule = *currentMolecule;
        
        // Calculate chain number
        int totNoSecs = tempMolecule->getSubsecSize(1);
        int chain_num = 1;
        while(change_section > (totNoSecs-1)) {
            chain_num++;
            totNoSecs = totNoSecs + tempMolecule->getSubsecSize(chain_num);
        }
        totNoSecs = totNoSecs - tempMolecule->getSubsecSize(chain_num) - 1;
        int local_section = change_section - totNoSecs - 1;
        
        // Try the change on temp molecule
        tempMolecule->changeMoleculeSingleMulti(local_section, chain_num);
        
        if (!tempMolecule->checkCalphas(chain_num)) {
            // Get coordinates from successful change
            auto coords = tempMolecule->getCoordinates();
            
            // Calculate total size needed
            int total_size = 0;
            for (auto& chain : coords) {
                total_size += chain.size() * 3;
            }
            
            // Allocate array
            double* data = new double[total_size];
            int idx = 0;
            
            // Copy coordinates to flat array
            for (auto& chain : coords) {
                for (auto& point : chain) {
                    data[idx++] = point.getX();
                    data[idx++] = point.getY();
                    data[idx++] = point.getZ();
                }
            }
            
            result.data = data;
            result.size = total_size;
            return result;
        }
    }
    
    // Only print message if we hit max attempts
    std::cout << "Failed to find valid conformation for section " << change_section 
            //   << " (local section " << change_section - totNoSecs - 1 
            //   << " in chain " << chain_num 
              << ") after " << MAX_ATTEMPTS << " attempts - CA distance check failed" << std::endl;
    
    return result;
}

    EXPORT bool accept_update() {
        if (!currentMolecule || !tempMolecule) return false;
        // Copy temp state to current
        *currentMolecule = *tempMolecule;
        return true;
    }


    EXPORT CoordArray get_current_coords() {
        CoordArray result = {nullptr, 0};
        if (!currentMolecule) return result;
        
        // Get coordinates
        auto coords = currentMolecule->getCoordinates();
        
        // Calculate total size needed
        int total_size = 0;
        for (auto& chain : coords) {
            total_size += chain.size() * 3;
        }
        
        // Allocate array
        double* data = new double[total_size];
        int idx = 0;
        
        // Copy coordinates to flat array
        for (auto& chain : coords) {
            for (auto& point : chain) {
                data[idx++] = point.getX();
                data[idx++] = point.getY();
                data[idx++] = point.getZ();
            }
        }
        
        result.data = data;
        result.size = total_size;
        return result;
    }

    EXPORT void free_coords(double* data) {
        if (data != nullptr) {
            delete[] data;
        }
    }

    EXPORT void cleanup_state() {
        if (currentMolecule) {
            delete currentMolecule;
            currentMolecule = nullptr;
        }
        if (tempMolecule) {
            delete tempMolecule;
            tempMolecule = nullptr;
        }
    }
}
