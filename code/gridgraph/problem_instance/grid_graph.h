
#ifndef TURNCOSTTOURSOPTIMIZATION_GRID_GRAPH_H
#define TURNCOSTTOURSOPTIMIZATION_GRID_GRAPH_H

#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/directed_graph.hpp>
#include "boost/multi_array.hpp"
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <stdlib.h>


typedef std::pair<int, int> coord;


namespace gg {

    class grid_graph {
    public:
        typedef boost::property<boost::vertex_name_t, coord> vertex_name_property_t;
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, vertex_name_property_t> boost_grid_graph_t;
        typedef boost_grid_graph_t::vertex_descriptor vertex_descriptor;

    private:
        //boost::property_map<boost_grid_graph_t, boost::vertex_name_t>::type vertex_name_map;
        std::map<coord, boost_grid_graph_t::vertex_descriptor> coord_to_vertex;

        void build_graph(std::set<coord> &coords);


    public:
        boost_grid_graph_t grid_graph_boost;

        inline boost_grid_graph_t::vertex_descriptor getVertex(coord c) const {
            return coord_to_vertex.at(c);
        }

        inline coord getCoord(boost_grid_graph_t::vertex_descriptor v) const {
            //vertex_name_map = boost::get(boost::vertex_name, grid_graph_boost);
            //assert(getVertex(boost::get(vertex_name_map, v))==v);
            return boost::get(boost::get(boost::vertex_name, grid_graph_boost), v);
        }

        grid_graph() { }

        grid_graph(std::set<coord> &coords);

        inline bool is_connected() const {
            std::vector<int> component(num_vertices(grid_graph_boost));
            int num = boost::connected_components(grid_graph_boost, &component[0]);
            return num == 1;
        }

        //Loads the grid graph from a file created by write_to_file
        void load_from_file(const std::string filename);

        void load_from_bitmap(const std::string filename, int threshold=100);

        //Saves the grid graph such that it can be reconstructed by load_from_file
        void write_to_file(const std::string filename);

        //bounding box
        int m_max_x;
        int m_min_x;
        int m_max_y;
        int m_min_y;

    };

}

#endif //TURNCOSTTOURSOPTIMIZATION_GRID_GRAPH_H
