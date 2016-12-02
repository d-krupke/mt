
#ifndef GEOM_ASPXOLVER_H
#define GEOM_APXSOLVER_H


#include <set>
#include <cstddef>
#include "../auxiliary.h"
#include "../cplex.hpp"
#include <boost/graph/properties.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

class ApproximateSolver
{
private:
    //instance
    const double m_DIST_COST;
    const double m_TURN_COST;
    const size_t m_NR_OF_ORIENTATIONS;
    const std::time_t m_DEADLINE;
    const std::set<coord> m_instance;

    //internal graph representation
    typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty> WeightedGraph;
    WeightedGraph instance_graph;
    using vd=WeightedGraph::vertex_descriptor;
    std::map<WeightedGraph::vertex_descriptor, conf> vertex_to_configuration;
    std::map<vd, vd> vertex_pair;
    std::map<coord, std::set<vd>> pi_confs;
    std::map<coord, coord> minimal_cycle_pair;

    //CPLEX
    IloEnv m_env;
    IloModel m_model;
    IloCplex m_cplex;

    std::map<WeightedGraph::edge_descriptor, IloNumVar> variables;

public:
    ApproximateSolver(const std::set<coord> instance, const size_t nr_of_orientations = 2, const double dist_weight = 0.0, const double turn_weight = 1.0, const unsigned long time_in_sec = 900);


    std::pair<double, std::vector<std::pair<coord, coord>>> solve(bool tour=true);
};


#endif //GEOM_ASPXOLVER_H
