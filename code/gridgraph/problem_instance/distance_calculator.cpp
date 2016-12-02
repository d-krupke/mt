
#include "distance_calculator.h"
gg::distance_calculator::distance_calculator(grid_graph &instance, double turn_cost, double dist_cost) 
	: TURN_COST{turn_cost}, DIST_COST{dist_cost}
{
        for (auto it = boost::vertices(instance.grid_graph_boost); it.first != it.second; ++it.first) {
            auto v = *it.first;
            north_map[v] = boost::add_vertex(dgraph);
            east_map[v] = boost::add_vertex(dgraph);
            south_map[v] = boost::add_vertex(dgraph);
            west_map[v] = boost::add_vertex(dgraph);
            boost::add_edge(north_map[v], east_map[v], TURN_COST, dgraph);
            boost::add_edge(east_map[v], north_map[v], TURN_COST, dgraph);
            boost::add_edge(east_map[v], south_map[v], TURN_COST, dgraph);
            boost::add_edge(south_map[v], east_map[v], TURN_COST, dgraph);
            boost::add_edge(south_map[v], west_map[v], TURN_COST, dgraph);
            boost::add_edge(west_map[v], south_map[v], TURN_COST, dgraph);
            boost::add_edge(west_map[v], north_map[v], TURN_COST, dgraph);
            boost::add_edge(north_map[v], west_map[v], TURN_COST, dgraph);
        }
        for (auto it = boost::edges(instance.grid_graph_boost); it.first != it.second; ++it.first) {
            auto e = *it.first;
            grid_graph::vertex_descriptor v = boost::source(e, instance.grid_graph_boost);
            grid_graph::vertex_descriptor w = boost::target(e, instance.grid_graph_boost);
            auto pos_v = instance.getCoord(v);
            auto pos_w = instance.getCoord(w);
            if (pos_v.first == pos_w.first) { //vertical
                if (pos_v.second < pos_w.second) { //v below w
                    boost::add_edge(north_map[v], north_map[w], DIST_COST, dgraph);
                    boost::add_edge(south_map[w], south_map[v], DIST_COST, dgraph);
                } else {
                    boost::add_edge(north_map[w], north_map[v], DIST_COST, dgraph);
                    boost::add_edge(south_map[v], south_map[w], DIST_COST, dgraph);
                }
            } else { //horizontal
                if (pos_v.first < pos_w.first) { //v left of w
                    boost::add_edge(east_map[v], east_map[w], DIST_COST, dgraph);
                    boost::add_edge(west_map[w], west_map[v], DIST_COST, dgraph);
                } else { //v right of w
                    boost::add_edge(west_map[v], west_map[w], DIST_COST, dgraph);
                    boost::add_edge(east_map[w], east_map[v], DIST_COST, dgraph);
                }
            }
        }
        auto num_vertices = boost::num_vertices(dgraph);
        distance_matrix.resize(boost::extents[num_vertices][num_vertices]);
        std::cout << "Calculate all distances..." << std::endl;
        boost::johnson_all_pairs_shortest_paths(dgraph, distance_matrix);
        std::cout << "Done." << std::endl;
    }


    double gg::distance_calculator::distance(grid_graph::vertex_descriptor v1, directions d_v1,
                                      grid_graph::vertex_descriptor v2,
                                      directions d_v2) {
        DGraph::vertex_descriptor aux_v1;
        switch (d_v1) {
            case directions::north:
                aux_v1 = north_map[v1];
                break;
            case directions::east:
                aux_v1 = east_map[v1];
                break;
            case directions::south:
                aux_v1 = south_map[v1];
                break;
            case directions::west:
                aux_v1 = west_map[v1];
                break;
            default: assert(false); break;
        }
        DGraph::vertex_descriptor aux_v2;
        switch (d_v2) {
            case directions::north:
                aux_v2 = south_map[v2];
                break;
            case directions::east:
                aux_v2 = west_map[v2];
                break;
            case directions::south:
                aux_v2 = north_map[v2];
                break;
            case directions::west:
                aux_v2 = east_map[v2];
                break;
            default: assert(false); break;
        }


        auto index_map = boost::get(boost::vertex_index, dgraph);
        auto v1_idx = boost::get(index_map, aux_v1);
        auto v2_idx = boost::get(index_map, aux_v2);
        return boost::get<double>(distance_matrix[v1_idx][v2_idx]);
    }
