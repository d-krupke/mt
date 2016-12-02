
#ifndef TURNCOSTTOURSOPTIMIZATION_INSTANCE_VISUALIZATION_H
#define TURNCOSTTOURSOPTIMIZATION_INSTANCE_VISUALIZATION_H

#include "grid_graph.h"

namespace gg {
    void visualize_instance(const std::string& filename, grid_graph& instance);

    void visualize_instance(const std::string& filename, grid_graph& instance, std::set<coord>& subset);

    void visualize_instance(const std::string& filename, grid_graph& instance, std::map<coord, double>& penalty);
}

#endif //TURNCOSTTOURSOPTIMIZATION_INSTANCE_VISUALIZATION_H
