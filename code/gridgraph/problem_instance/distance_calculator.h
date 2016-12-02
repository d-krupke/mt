
#ifndef TURNCOSTTOURSOPTIMIZATION_DISTANCE_CALCULATOR_H
#define TURNCOSTTOURSOPTIMIZATION_DISTANCE_CALCULATOR_H

#include <boost/multi_array.hpp>
#include <boost/graph/directed_graph.hpp>
#include "grid_graph.h"

namespace gg{

    enum class directions {
        north, east, south, west
    };

    class distance_calculator {

    private:
        boost::multi_array<double, 2> distance_matrix;


        //integral edge weights
        typedef boost::property<boost::edge_weight_t, double> FloatEdgeWeightProperty;
        //Link vertices to their origin field and their direction
        typedef boost::property<boost::vertex_name_t, grid_graph::vertex_descriptor> MapBackVertexProperty;
        //Directed Graph with the two above properties
        typedef boost::directed_graph<MapBackVertexProperty, FloatEdgeWeightProperty> DGraph;

        //std::map<DGraph::vertex_descriptor, std::map<DGraph::vertex_descriptor, int>> distance_matrix;

        std::map<grid_graph::vertex_descriptor, DGraph::vertex_descriptor> north_map;
        std::map<grid_graph::vertex_descriptor, DGraph::vertex_descriptor> east_map;
        std::map<grid_graph::vertex_descriptor, DGraph::vertex_descriptor> south_map;
        std::map<grid_graph::vertex_descriptor, DGraph::vertex_descriptor> west_map;

        DGraph dgraph;

		const double TURN_COST;
		const double DIST_COST;

    public:


        distance_calculator(grid_graph &instance, double turn_cost=1, double dist_cost=0);


        double distance(grid_graph::vertex_descriptor v1, directions d_v1, grid_graph::vertex_descriptor v2, directions d_v2);
    };

}
#endif //TURNCOSTTOURSOPTIMIZATION_DISTANCE_CALCULATOR_H
