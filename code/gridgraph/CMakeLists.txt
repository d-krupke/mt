cmake_minimum_required(VERSION 3.6)
project(gridgraph)


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


add_executable(fc     fc/base_solver.h
        fc/formulation_1.cpp
        fc/formulation_1.h
        fc/formulation_2.cpp
        fc/formulation_2.h
        fc/formulation_3.cpp
        fc/formulation_3.h
        fc/main.cpp
        fc/solution_printer.cpp
        fc/solution_printer.h
        problem_instance/distance_calculator.cpp
        problem_instance/distance_calculator.h
        problem_instance/grid_graph.cpp
        problem_instance/grid_graph.h
        problem_instance/save_and_load.cpp
        problem_instance/save_and_load.h)
add_executable(subset       subset/base_solver.h
        subset/formulation_1.cpp
        subset/formulation_1.h
        subset/formulation_2.cpp
        subset/formulation_2.h
        subset/formulation_3.cpp
        subset/formulation_3.h
        subset/main.cpp
        problem_instance/distance_calculator.cpp
        problem_instance/distance_calculator.h
        problem_instance/grid_graph.cpp
        problem_instance/grid_graph.h
        problem_instance/instance_visualization.cpp
        problem_instance/instance_visualization.h
        problem_instance/save_and_load.cpp
        problem_instance/save_and_load.h
)
add_executable(penalty      penalty/formulation_1.cpp
        penalty/formulation_1.h
        penalty/formulation_2.cpp
        penalty/formulation_2.h
        penalty/formulation_3.cpp
        penalty/formulation_3.h
        penalty/main.cpp
        problem_instance/distance_calculator.cpp
        problem_instance/distance_calculator.h
        problem_instance/grid_graph.cpp
        problem_instance/grid_graph.h
        problem_instance/instance_visualization.cpp
        problem_instance/instance_visualization.h
        problem_instance/save_and_load.cpp
        problem_instance/save_and_load.h
        )
add_executable(instance_creator
        ./problem_instance/main.cpp
        ./problem_instance/grid_graph.cpp
        ./problem_instance/save_and_load.cpp
        ./problem_instance/random_instances.cpp)
target_include_directories(fc   PRIVATE ${INCLUDE_DIRS})
target_link_libraries(fc    ${LIBRARIES})
target_include_directories(subset  PRIVATE ${INCLUDE_DIRS})
target_link_libraries(subset    ${LIBRARIES})
target_include_directories(penalty  PRIVATE ${INCLUDE_DIRS})
target_link_libraries(penalty    ${LIBRARIES})
