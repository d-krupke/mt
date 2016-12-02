cmake_minimum_required(VERSION 3.6)
project(geom)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
# Boost and its components
find_package( Boost COMPONENTS program_options graph REQUIRED )
find_package( CPLEX REQUIRED )
find_package( Threads REQUIRED)
if ( NOT Boost_FOUND )
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(LIBRARIES ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} Threads::Threads)
set(INCLUDE_DIRS ${CPLEX_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})


set(SOURCE_FILES main.cpp auxiliary.h ip/OptimalSolver.cpp ip/OptimalSolver.h visualization/Visualization.cpp visualization/Visualization.h aa/ApproximateSolver.h aa/ApproximateSolver.cpp)
add_executable(geom ${SOURCE_FILES})

target_include_directories(geom  PRIVATE ${INCLUDE_DIRS})
target_link_libraries(geom   ${LIBRARIES})
