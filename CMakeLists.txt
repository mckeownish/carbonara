cmake_minimum_required(VERSION 3.10)
project(Carbonara)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Force use of Conda Python
set(Python3_ROOT_DIR "/opt/anaconda3/envs/NLE_lab")
set(Python3_EXECUTABLE "/opt/anaconda3/envs/NLE_lab/bin/python")

# Find Python
find_package(Python3 COMPONENTS Development REQUIRED)

# Print Python info
message(STATUS "Found Python at: ${Python3_EXECUTABLE}")
message(STATUS "Using Python include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "Using Python libraries: ${Python3_LIBRARIES}")

# Include directories
include_directories(src)
include_directories(${Python3_INCLUDE_DIRS})

# Source files
set(SOURCES
    src/point.cpp
    src/polyHelix.cpp
    src/randomMolGen.cpp
    src/ktlMoleculeRandom.cpp
    src/experimentalData.cpp
    src/hydrationShellRandom.cpp
    src/skmt.cpp
    src/writheFP.cpp
    src/moleculeFitAndState.cpp
    src/Logger.cpp
    src/helpers.cpp
    src/parameters.cpp
)

# Main executables
add_executable(predictStructureQvary ${SOURCES} src/mainPredictionFinalQvar.cpp)
add_executable(generate_structure ${SOURCES} src/Flexible_generator.cpp)
add_executable(single_fit ${SOURCES} src/single_fit.cpp)
add_executable(AppWAXSiS ${SOURCES} src/newScatteringTestBench.cpp)


# Python interface library
add_library(molecule_interface SHARED
    ${SOURCES}
    src/molecule_interface.cpp
)

# Set properties for the interface library
set_target_properties(molecule_interface PROPERTIES 
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    C_VISIBILITY_PRESET default
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
)

if(APPLE)
    target_compile_options(molecule_interface PRIVATE -fvisibility=default)
    set_target_properties(molecule_interface PROPERTIES 
        LINK_FLAGS "-undefined dynamic_lookup"
        MACOSX_RPATH ON
    )
endif()

# Link Python to the interface library
target_link_libraries(molecule_interface ${Python3_LIBRARIES})

# Set output directory for executables
set_target_properties(predictStructureQvary generate_structure single_fit
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)