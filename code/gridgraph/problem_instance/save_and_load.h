
#ifndef TURNCOSTTOURSOPTIMIZATION_SAVE_AND_LOAD_H
#define TURNCOSTTOURSOPTIMIZATION_SAVE_AND_LOAD_H


#include "grid_graph.h"

namespace gg{
    void save_instance(grid_graph& instance, std::string filename);

    grid_graph load_instance(std::string filename);

    void save_instance(grid_graph& instance, std::set<coord>& subset, std::string filename);

    grid_graph load_instance(std::string filename, std::set<coord>* subset);

    void save_instance(grid_graph& instance, std::map<coord, double>& penalty, std::string filename);

    grid_graph load_instance(std::string filename, std::map<coord, double>& penalty);
}


#endif //TURNCOSTTOURSOPTIMIZATION_SAVE_AND_LOAD_H
