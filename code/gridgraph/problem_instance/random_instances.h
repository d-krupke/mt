
#ifndef TURNCOSTTOURSOPTIMIZATION_RANDOM_INSTANCES_H
#define TURNCOSTTOURSOPTIMIZATION_RANDOM_INSTANCES_H

#include <stdlib.h>
#include <random>
#include "grid_graph.h"
//#include "png++/png.hpp"

namespace gg{
    grid_graph create_random_grid_graph_from_png(std::string filename, size_t desired_size, float thinout_fac, float inc_direct_nbr, float inc_diag_nbr);

    std::set<coord> create_random_subset_from_png(grid_graph& instance, std::string filename, int desired_subset_size, float margin);

    std::map<coord, double> create_random_penalty_from_png(grid_graph& instance, std::string filename, double mean=2, double stddev=1);

    std::set<coord> create_random_subset(grid_graph& instance, size_t desired_size);

    std::map<coord, double> create_random_penalty(grid_graph& instance, double mean=2, double stddev=1);

    grid_graph create_random_grid_graph(size_t desired_size);
    grid_graph create_random_grid_graph(size_t desired_size, int max_x, int max_y, float delete_prob_perc=0.01f, float inc_direct_nbr=2.0f, float inc_diag_nbr=1.0f);
}

#endif //TURNCOSTTOURSOPTIMIZATION_RANDOM_INSTANCES_H
