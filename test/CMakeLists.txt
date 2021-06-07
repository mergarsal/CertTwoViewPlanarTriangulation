project(projection_relpose-Tests CXX)

# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})

# required by SDPA
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")



# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})



add_library(exp_help_test 
                experimentsHelper.h 
                ../utils/triangulationPlanar.h
                ../utils/generatePointCloudPlanar.h
                #                
                experimentsHelper.cpp 
                ../utils/triangulationPlanar.cpp
                ../utils/generatePointCloudPlanar.cpp)





## Generic test synthetic
add_executable(generic_test ${CMAKE_CURRENT_SOURCE_DIR}/generic_test.cpp)
target_link_libraries(generic_test  
        OptTwoPlanarTriangulation  
        exp_help_test

)
        



