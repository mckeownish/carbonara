cmake_minimum_required(VERSION 3.10)
project(MinimalTest)

# Force use of Conda Python
set(Python3_ROOT_DIR "/opt/anaconda3/envs/NLE_lab")
set(Python3_EXECUTABLE "/opt/anaconda3/envs/NLE_lab/bin/python")

find_package(Python3 COMPONENTS Development REQUIRED)

add_library(minimal_test SHARED src/minimal_test.cpp)

# Explicitly set compiler flags for macOS
if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
    set_target_properties(minimal_test PROPERTIES
        LINK_FLAGS "-undefined dynamic_lookup"
        MACOSX_RPATH ON
    )
endif()

target_include_directories(minimal_test PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(minimal_test ${Python3_LIBRARIES})

set_target_properties(minimal_test PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)