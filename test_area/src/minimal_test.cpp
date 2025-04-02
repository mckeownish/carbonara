#include <iostream>

// Explicit symbol export for macOS
#ifdef __APPLE__
    #define EXPORT __attribute__((visibility("default")))
#else
    #define EXPORT
#endif

extern "C" {
    // Simple array return structure
    struct CoordArray {
        double* data;
        int size;
    };

    EXPORT CoordArray create_test_coords() {
        std::cout << "Creating test coords" << std::endl;
        int size = 9;  // 3 points x 3 coordinates
        double* coords = new double[size];
        
        // Initialize array manually
        coords[0] = 1.0; coords[1] = 2.0; coords[2] = 3.0;  // Point 1
        coords[3] = 4.0; coords[4] = 5.0; coords[5] = 6.0;  // Point 2
        coords[6] = 7.0; coords[7] = 8.0; coords[8] = 9.0;  // Point 3
        
        CoordArray result;
        result.data = coords;
        result.size = size;
        
        std::cout << "Created array of size " << size << std::endl;
        return result;
    }

    EXPORT void free_coords(double* data) {
        delete[] data;
    }
}