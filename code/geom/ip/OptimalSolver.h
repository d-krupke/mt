
#ifndef GEOM_OPTIMALSOLVER_H
#define GEOM_OPTIMALSOLVER_H


#include <set>
#include <cstddef>
#include "../auxiliary.h"
#include "../cplex.hpp"
#include <boost/graph/properties.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/connected_components.hpp>

class OptimalSolver
{
private:
    //instance
    const double m_DIST_COST;
    const double m_TURN_COST;
    const std::time_t m_DEADLINE;
    const std::set<coord> m_instance;

    //internal graph representation
    typedef boost::property<boost::edge_weight_t, double> internal_graph_edge_weight_property;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, internal_graph_edge_weight_property> internal_graph;
    internal_graph graph;
    using vd=internal_graph::vertex_descriptor;
    std::map<internal_graph::vertex_descriptor, conf> vertex_to_configuration;
    std::map<vd, vd> vertex_pair;
    std::map<coord, std::set<vd>> pi_confs;
    std::map<coord, coord> minimal_cycle_pair;

    //CPLEX
    IloEnv m_env;
    IloModel m_model;
    IloCplex m_cplex;

    std::map<internal_graph::edge_descriptor, IloBoolVar> variables;
    bool cplex_params_set = false;

public:
    OptimalSolver(const std::set<coord> instance, const size_t nr_of_orientations = 2, bool cc = false, const double dist_weight = 0.0, const double turn_weight = 1.0, const std::time_t time_in_sec = 900) :
            m_DIST_COST{dist_weight}, m_TURN_COST{turn_weight}, m_DEADLINE{std::time(nullptr) + time_in_sec}, m_instance{instance}, m_model{m_env}
    {
        std::cout << "Nr of orientations="<<nr_of_orientations<<" cc="<<cc<<" dist_weight="<<dist_weight<<" turn_weight="<<turn_weight<<std::endl;
        auto angle_steps = M_PI / nr_of_orientations;

        //create vertices
        for (auto p: instance) {
            for (auto i = 0; i < nr_of_orientations; ++i) {
                auto heading = i * angle_steps;
                auto counter_heading = heading + M_PI;

                conf a{p, heading};
                conf b{p, counter_heading};
                auto v1 = boost::add_vertex(graph);
                vertex_to_configuration[v1] = a;
                auto v2 = boost::add_vertex(graph);
                vertex_to_configuration[v2] = b;
                vertex_pair[v1] = v2;
                vertex_pair[v2] = v1;
                pi_confs[p].insert(v1);
            }
        }

        IloExpr obj_expr(m_env);

        //create edges
        for (auto v_it = boost::vertices(graph); v_it.first != v_it.second; ++v_it.first) {
            auto v = *v_it.first;
            for (auto w_it = boost::vertices(graph); w_it.first != v_it.second; ++w_it.first) {
                auto w = *w_it.first;
                if (vertex_to_configuration[v].first == vertex_to_configuration[w].first) continue;
                if (v < w) continue; // only vw or wv, not both

                auto cost = costs(vertex_to_configuration[v], vertex_to_configuration[w], m_DIST_COST, m_TURN_COST);
                auto e = boost::add_edge(v, w, cost, graph);
                assert(e.second);

                //add edge to cplex
                auto var = IloBoolVar(m_env);
                variables[e.first] = var;
                obj_expr += cost * var;
            }
        }

        if (cc) {
            //Add minimal cycle
            for (auto p: instance) {
                double closest;
                bool closest_initialized = false;
                for (auto q: instance) {
                    if (p == q) continue;
                    auto distc = dist(p, q);
                    if (closest > distc || !closest_initialized) {
                        closest = distc;
                        closest_initialized = true;
                        minimal_cycle_pair[p] = q;
                    }
                }
                assert(closest_initialized);
                auto v = *pi_confs[p].begin();
                auto w = vertex_pair[v];
                double cost = closest * 2 * m_DIST_COST + M_PI * 2 * m_TURN_COST;
                auto e = boost::add_edge(v, w, cost, graph).first;
                variables[e] = IloBoolVar(m_env);
                obj_expr+=cost*variables[e];
            }
        }

        IloObjective objective = IloMinimize(m_env, obj_expr);
        m_model.add(objective);

        //Add constraints
        for (auto p: instance) {
            IloExpr visited_once{m_env};

            for (auto v: pi_confs[p]) {
                auto w = vertex_pair[v];
                IloExpr equal_vists{m_env};

                for (auto e_it = boost::out_edges(v, graph); e_it.first != e_it.second; ++e_it.first) {
                    auto e = *e_it.first;
                    assert(variables.count(e) > 0);
                    auto var = variables[e];
                    visited_once += 1 * var;
                    equal_vists += 1 * var;
                }
                for (auto e_it = boost::out_edges(w, graph); e_it.first != e_it.second; ++e_it.first) {
                    auto e = *e_it.first;
                    auto var = variables[e];
                    assert(variables.count(e) > 0);
                    equal_vists += (-1) * var;
                }

                m_model.add(IloRange{m_env, 0, equal_vists, 0});

            }

            m_model.add(IloRange{m_env, 1, visited_once, 1});
        }

        m_cplex = IloCplex{m_model};
        //std::cout << m_model <<std::endl;
    }

    size_t separate_subcycles()
    {
        typedef boost::property<boost::vertex_name_t, coord> map_back_property;
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, map_back_property, boost::no_property> comp_graph;
        std::map<coord, comp_graph::vertex_descriptor> gg_to_cg;

        comp_graph cg;
        size_t cuts_added = 0;

        //add vertices
        auto map_back = boost::get(boost::vertex_name, cg);
        for (auto p: m_instance) {
            auto v_cg = boost::add_vertex(cg);
            boost::put(map_back, v_cg, p);
            gg_to_cg[p] = v_cg;
        }

        //add edges
        std::map<internal_graph::edge_descriptor, IloNum> solution;
        for (auto ev: variables) {
            solution[ev.first] = m_cplex.getValue(ev.second);
        }
        for (auto val_it: solution) {
            if (val_it.second > 0.1) {
                auto v_gg = vertex_to_configuration[boost::target(val_it.first, graph)].first;
                auto w_gg = vertex_to_configuration[boost::source(val_it.first, graph)].first;
                if (v_gg < w_gg) {
                    boost::add_edge(gg_to_cg[v_gg], gg_to_cg[w_gg], cg);
                }
            }
        }

        //calculate components
        std::vector<int> component(boost::num_vertices(cg));
        int num = boost::connected_components(cg, &component[0]);
        if (num > 1) {
            //create map with all components
            std::map<int, std::vector<coord>> comps;
            for (auto v_it = boost::vertices(cg); v_it.first != v_it.second; ++v_it.first) {
                auto v = *v_it.first;
                auto cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v)];
                assert(cmp_idx <= num);
                if (!comps.count(cmp_idx)) comps[cmp_idx] = std::vector<coord>();
                comps[cmp_idx].push_back(boost::get(boost::get(boost::vertex_name, cg), v));
            }

            //add constraint for each component
            for (auto comp: comps) {
                std::cout << "Comp " << comp.first << " with " << comp.second.size() << " elements" << std::endl;
                IloExpr constr{m_env};
                for (auto p: comp.second) {
                    for (auto v: pi_confs[p]) {
                        auto w = vertex_pair[v];
                        for (auto e_it = boost::out_edges(v, graph); e_it.first != e_it.second; ++e_it.first) {
                            auto e = *e_it.first;
                            //Check if it connects two comps
                            auto v1 = vertex_to_configuration[boost::target(e, graph)].first;
                            auto v2 = vertex_to_configuration[boost::source(e, graph)].first;
                            auto v1_aux = gg_to_cg[v1];
                            auto v2_aux = gg_to_cg[v2];
                            auto v1_cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v1_aux)];
                            auto v2_cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v2_aux)];
                            if (v1_cmp_idx != v2_cmp_idx) {
                                assert(solution[*e_it.first] < 0.1);
                                constr += 1 * variables[*e_it.first];
                            }
                        }
                        for (auto e_it = boost::out_edges(w, graph); e_it.first != e_it.second; ++e_it.first) {
                            auto e = *e_it.first;
                            //check if it connects two comps
                            auto v1 = vertex_to_configuration[boost::target(e, graph)].first;
                            auto v2 = vertex_to_configuration[boost::source(e, graph)].first;
                            auto v1_aux = gg_to_cg[v1];
                            auto v2_aux = gg_to_cg[v2];
                            auto v1_cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v1_aux)];
                            auto v2_cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v2_aux)];
                            if (v1_cmp_idx != v2_cmp_idx) {
                                assert(solution[*e_it.first] < 0.1);
                                constr += 1 * variables[*e_it.first];
                            }
                        }
                    }
                }
                m_model.add(IloRange(m_env, 2, constr, IloInfinity));
                ++cuts_added;
            }
            std::cout << "Solution consists of " << num << " components. Added " << cuts_added << " cuts" << std::endl;
            assert(cuts_added > 0);
        } else {
            std::cout << "Solution is connected. No cuts added." << std::endl;
        }
        return cuts_added;
    }

    IloAlgorithm::Status solve()
    {
        if (!cplex_params_set) {
            // m_cplex.setParam(IloCplex::NumericalEmphasis, 1);
            m_cplex.setParam(IloCplex::EpInt, 0);
            m_cplex.setParam(IloCplex::EpGap, 0);
            m_cplex.setParam(IloCplex::EpOpt, 1e-9);
            m_cplex.setParam(IloCplex::EpAGap, 0);
            cplex_params_set = true;
        }
        if (m_DEADLINE != 0) {
            std::time_t now = std::time(nullptr);
            long resttime = m_DEADLINE - now;
            if (resttime <= 0) {
                std::cout << "Abort solve because no time left." << std::endl;
                throw (-1);
            }
            m_cplex.setParam(IloCplex::TiLim, resttime);
        }
        std::cout << "Solving IP...." << std::endl;
        //Solve
        try {
            if (!m_cplex.solve()) {
                m_env.error() << "Failed to optimize LP: " << m_cplex.getStatus() << std::endl;
                throw (-1);
            }
        } catch (IloException e) {
            std::cerr << "Error while solving: " << e << std::endl;
        }
        std::cout << "Solved IP to optimality with solution: " << m_cplex.getObjValue() << std::endl;
        if (m_DEADLINE < std::time(nullptr)) throw (-1);
        return m_cplex.getStatus();
    }

    double get_objective_value()
    {
        constexpr double RETURN_VALUE_IF_OUT_OF_TIME = -1;
        if (m_DEADLINE < std::time(nullptr)) return RETURN_VALUE_IF_OUT_OF_TIME;
        double ret;
        try {
            ret = m_cplex.getObjValue();
        } catch (...) {
            std::cerr << "ERROR: Exception thrown while querying objective value. Return 0." << std::endl;
            return RETURN_VALUE_IF_OUT_OF_TIME;
        }
        return ret;
    }

    std::vector<std::pair<coord, coord>> get_solution()
    {
        std::vector<std::pair<coord, coord>> solution;
        for (auto e_it = boost::edges(graph); e_it.first != e_it.second; ++e_it.first) {
            auto e = *e_it.first;
            auto v = boost::source(e, graph);
            auto w = boost::target(e, graph);
            if (m_cplex.getValue(variables[e]) > 0.1) {
                solution.push_back({vertex_to_configuration[v].first, vertex_to_configuration[w].first});
            }
        }
        return solution;
    }

    std::vector<coord> get_tour()
    {
        assert(false);//TODO
    }
};


#endif //GEOM_OPTIMALSOLVER_H
