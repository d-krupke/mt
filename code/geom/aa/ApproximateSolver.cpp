
#include "ApproximateSolver.h"

ApproximateSolver::ApproximateSolver(const std::set<coord> instance, const size_t nr_of_orientations, const double dist_weight, const double turn_weight, const unsigned long time_in_sec) :
        m_DIST_COST{dist_weight}, m_TURN_COST{turn_weight}, m_NR_OF_ORIENTATIONS{nr_of_orientations}, m_DEADLINE{std::time(nullptr) + time_in_sec}, m_instance{instance}, m_model{m_env}
{
    auto angle_steps = M_PI / nr_of_orientations;

    std::cout << "Creating vertices..." << std::endl;
    //create vertices
    for (auto p: instance) {
        for (auto i = 0; i < nr_of_orientations; ++i) {
            auto heading = i * angle_steps;
            auto counter_heading = heading + M_PI;

            conf a{p, heading};
            conf b{p, counter_heading};
            auto v1 = boost::add_vertex(instance_graph);
            vertex_to_configuration[v1] = a;
            auto v2 = boost::add_vertex(instance_graph);
            vertex_to_configuration[v2] = b;
            vertex_pair[v1] = v2;
            vertex_pair[v2] = v1;
            pi_confs[p].insert(v1);
        }
    }

    std::cout << "Creating Edges..." << std::endl;
    IloExpr obj_expr(m_env);
    //create edges
    for (auto v_it = boost::vertices(instance_graph); v_it.first != v_it.second; ++v_it.first) {
        auto v = *v_it.first;
        for (auto w_it = boost::vertices(instance_graph); w_it.first != v_it.second; ++w_it.first) {
            auto w = *w_it.first;
            if (vertex_to_configuration[v].first == vertex_to_configuration[w].first) continue;
            if (v < w) continue; // only vw or wv, not both

            auto cost = costs(vertex_to_configuration[v], vertex_to_configuration[w], m_DIST_COST, m_TURN_COST);
            auto e = boost::add_edge(v, w, cost, instance_graph);
            assert(e.second);

            //add edge to cplex
            auto var = IloNumVar(m_env, 0.0, 1.0);
            variables[e.first] = var;
            obj_expr += cost * var;
        }
    }

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
        for(auto v: pi_confs[p]){
            auto w = vertex_pair[v];
            double cost = closest * 2 * m_DIST_COST + M_PI * 2 * m_TURN_COST;
            auto e = boost::add_edge(v, w, cost, instance_graph).first;
            variables[e] = IloNumVar(m_env, 0.0, 1.0);
            obj_expr+=cost*variables[e];
        }
    }

    IloObjective objective = IloMinimize(m_env, obj_expr);
    m_model.add(objective);

    std::cout << "Adding constraints..." << std::endl;
    //Add constraints
    for (auto p: instance) {
        IloExpr visited_once{m_env};

        for (auto v: pi_confs[p]) {
            auto w = vertex_pair[v];
            IloExpr equal_vists{m_env};

            for (auto e_it = boost::out_edges(v, instance_graph); e_it.first != e_it.second; ++e_it.first) {
                auto e = *e_it.first;
                assert(variables.count(e) > 0);
                auto var = variables[e];
                visited_once += 1 * var;
                equal_vists += 1 * var;
            }
            for (auto e_it = boost::out_edges(w, instance_graph); e_it.first != e_it.second; ++e_it.first) {
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

std::pair<double, std::vector<std::pair<coord, coord>>> ApproximateSolver::solve(bool tour)
{

    m_cplex.setParam(IloCplex::EpInt, 0);
    m_cplex.setParam(IloCplex::EpGap, 0);
    m_cplex.setParam(IloCplex::EpOpt, 1e-9);
    m_cplex.setParam(IloCplex::EpAGap, 0);

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
    //**Build matching graph**
    //Get vertices
    std::set<vd> selected_vertices;
    for (auto p: m_instance) {
        bool inited = false;
        double best_value = 0.0;
        vd best_vertex;
        for (auto v: pi_confs[p]) {
            auto sum_edge = 0.0;
            for (auto e_it = boost::out_edges(v, instance_graph); e_it.first != e_it.second; ++e_it.first) {
                auto e = *e_it.first;
                auto val = m_cplex.getValue(variables[e]);
                sum_edge += val;
            }
            if(!inited || sum_edge>best_value){
                inited = true;
                best_value = sum_edge;
                best_vertex = v;
            }
        }
        selected_vertices.insert(best_vertex);
        selected_vertices.insert(vertex_pair[best_vertex]);
        assert(inited);
    }

    //Create Graph
    WeightedGraph matching_graph;
    std::map<vd, vd> instance_to_matching_graph;
    std::map<vd, vd> matching_to_instance_graph;
    for (auto v: selected_vertices) {
        auto v_sub = boost::add_vertex(matching_graph);
        instance_to_matching_graph[v] = v_sub;
        matching_to_instance_graph[v_sub] = v;
    }
    for (auto v: selected_vertices) {
        for (auto w: selected_vertices) {
            if (v <= w) continue;
            auto e = boost::edge(v, w, instance_graph);
            if (!e.second) {
                continue;
            }
            auto c = boost::get(boost::get(boost::edge_weight, instance_graph), e.first);
            auto e_sub = boost::add_edge(instance_to_matching_graph[v], instance_to_matching_graph[w], c, matching_graph);
        }
    }

    //**Solve MinPerfMatching**
    IloEnv matching_environment;
    IloModel matching_model{matching_environment};
    std::map<WeightedGraph::edge_descriptor, IloBoolVar> matching_variables;
    IloExpr obj_expr{matching_environment};
    for(auto e_it=boost::edges(matching_graph); e_it.first!=e_it.second; ++e_it.first){
        auto var = IloBoolVar{matching_environment};
        auto e = *e_it.first;
        matching_variables[e] = var;
        obj_expr+= boost::get(boost::get(boost::edge_weight, matching_graph), e)*var;
    }
    matching_model.add(IloMinimize(matching_environment, obj_expr));
    for(auto v=boost::vertices(matching_graph); v.first!=v.second; ++v.first){
        IloExpr constr{matching_environment};
        assert(boost::degree(*v.first, matching_graph)>0);
        for(auto e_it = boost::out_edges(*v.first, matching_graph); e_it.first!=e_it.second; ++e_it.first){
            constr+=matching_variables[*e_it.first];
        }
        matching_model.add(IloRange{matching_environment, 1, constr, 1});
    }
    IloCplex matching_cplex{matching_model};
    //std::cout << matching_model;
    matching_cplex.solve();

    //Solution Edges
    std::set<WeightedGraph::edge_descriptor> solution_edges_in_matching_graph;
    for(auto e_it=boost::edges(matching_graph); e_it.first!=e_it.second; ++e_it.first) {
        auto var = matching_variables[*e_it.first];
        if(matching_cplex.getValue(var)<0.1) continue;
        solution_edges_in_matching_graph.insert(*e_it.first);
    }

    //Create Solution Graph
    WeightedGraph solution_graph;
    std::map<coord, vd> coord_to_solution_vertex;
    std::map<vd, coord> solution_vertex_to_coord;
    std::map<WeightedGraph::edge_descriptor, WeightedGraph::edge_descriptor> solution_edge_to_matching_edge;
    for(auto p: m_instance) {
        auto v = boost::add_vertex(solution_graph);
        coord_to_solution_vertex[p] = v;
        solution_vertex_to_coord[v] = p;
    }
    for(auto e_it=boost::edges(matching_graph); e_it.first!=e_it.second; ++e_it.first) {
        auto var = matching_variables[*e_it.first];
        if(matching_cplex.getValue(var)<0.1) continue;
        auto conf_v = vertex_to_configuration[matching_to_instance_graph[boost::source(*e_it.first, matching_graph)]];
        auto v = coord_to_solution_vertex[conf_v.first];
        auto conf_w = vertex_to_configuration[matching_to_instance_graph[boost::target(*e_it.first, matching_graph)]];
        auto w = coord_to_solution_vertex[conf_w.first];
        auto e = boost::add_edge(v, w, boost::get(boost::get(boost::edge_weight, matching_graph), *e_it.first), solution_graph).first;
        solution_edge_to_matching_edge[e] = *e_it.first;
    }

    //If cycle cover demanded, extract solution and return
    if(!tour){
        double cost = 0.0;
        std::vector<std::pair<coord, coord>> solution;
        for(auto e_it=boost::edges(solution_graph); e_it.first!=e_it.second; ++e_it.first) {
            solution.push_back({solution_vertex_to_coord[boost::source(*e_it.first, solution_graph)], solution_vertex_to_coord[boost::target(*e_it.first, solution_graph)]});
            cost += boost::get(boost::get(boost::edge_weight, solution_graph), *e_it.first);
        }
        return {cost, solution};
    }

    //Calculate Components in Solution Graph
    std::vector<int> component(boost::num_vertices(solution_graph));
    int num = boost::connected_components(solution_graph, &component[0]);
    std::map<int, std::vector<vd>> components;
    for(auto v_it=boost::vertices(solution_graph); v_it.first!=v_it.second; ++v_it.first){
        auto v = *v_it.first;
        auto cmp_idx = component[boost::get(boost::get(boost::vertex_index, solution_graph), v)];
        components[cmp_idx].push_back(v);
    }

    //Component Graph to calculate MST
    WeightedGraph mst_graph;
    std::map<int, vd> cmp_id_to_vertex;
    std::map<vd, int> vertex_to_cmp_id;
    for(auto comp: components){
        auto v = boost::add_vertex(mst_graph);
        cmp_id_to_vertex[comp.first] = v;
        vertex_to_cmp_id[v] = comp.first;
    }
    std::map<WeightedGraph::edge_descriptor, std::pair<coord, coord>> mst_edge_to_coords;
    for(auto comp_a: components) {
        for(auto comp_b: components) {
            if(comp_a.first<=comp_b.first) continue;
            //Calculate cheapest connection
            std::pair<coord, coord> best_connection{solution_vertex_to_coord[*comp_a.second.begin()], solution_vertex_to_coord[*comp_b.second.begin()]};
            double best_cost = dist(best_connection.first, best_connection.second);
            for(auto p: comp_a.second) {
                for(auto q: comp_b.second) {
                    auto cost = dist(solution_vertex_to_coord[p], solution_vertex_to_coord[q]);
                    if(cost<best_cost) {
                        best_connection = {solution_vertex_to_coord[p], solution_vertex_to_coord[q]};
                        best_cost = cost;
                    }
                }
            }

            auto e = boost::add_edge(cmp_id_to_vertex[comp_a.first], cmp_id_to_vertex[comp_b.first], best_cost, mst_graph);
            mst_edge_to_coords[e.first] = best_connection;
        }
    }
    std::vector < WeightedGraph::edge_descriptor > spanning_tree;
    boost::kruskal_minimum_spanning_tree(mst_graph, std::back_inserter(spanning_tree));

    //Connect Components
    for(auto e: spanning_tree){
        auto coords = mst_edge_to_coords[e];
        auto solution_v = coord_to_solution_vertex[coords.first];
        auto solution_w = coord_to_solution_vertex[coords.second];

        std::vector<vd> v_vertices;
        for(auto e_it=boost::out_edges(solution_v, solution_graph); e_it.first!=e_it.second; ++e_it.first){
            auto e_m = solution_edge_to_matching_edge[*e_it.first];
            auto v_a = boost::target(e_m, matching_graph);
            auto v_b = boost::source(e_m, matching_graph);
            auto coord_a = vertex_to_configuration[matching_to_instance_graph[v_a]].first;
            auto coord_b = vertex_to_configuration[matching_to_instance_graph[v_b]].first;
            if(coord_a == coords.first) v_vertices.push_back(v_a);
            if(coord_b == coords.first) v_vertices.push_back(v_b);
        }
        std::vector<vd> w_vertices;
        for(auto e_it=boost::out_edges(solution_w, solution_graph); e_it.first!=e_it.second; ++e_it.first){
            auto e_m = solution_edge_to_matching_edge[*e_it.first];
            auto v_a = boost::target(e_m, matching_graph);
            auto v_b = boost::source(e_m, matching_graph);
            auto coord_a = vertex_to_configuration[matching_to_instance_graph[v_a]].first;
            auto coord_b = vertex_to_configuration[matching_to_instance_graph[v_b]].first;
            if(coord_a == coords.second) w_vertices.push_back(v_a);
            if(coord_b == coords.second) w_vertices.push_back(v_b);
        }
        assert(v_vertices.size()==2 && w_vertices.size()==2);


        auto cost_map = boost::get(boost::edge_weight, matching_graph);
        auto e1 = boost::edge(v_vertices[0], w_vertices[0], matching_graph).first;
        auto e2 = boost::edge(v_vertices[1], w_vertices[1], matching_graph).first;
        std::pair<WeightedGraph::edge_descriptor, WeightedGraph::edge_descriptor> best_pair{e1, e2};
        double best_pair_value = boost::get(cost_map, e1)+boost::get(cost_map, e2);
        auto e3 = boost::edge(v_vertices[0], w_vertices[1], matching_graph).first;
        auto e4 = boost::edge(v_vertices[0], w_vertices[1], matching_graph).first;
        if(boost::get(cost_map, e3)+boost::get(cost_map, e4)<best_pair_value) best_pair = {e3, e4};
        solution_edges_in_matching_graph.insert(best_pair.first);
        solution_edges_in_matching_graph.insert(best_pair.second);
    }

    double cost = 0.0;
    std::vector<std::pair<coord, coord>> solution;
    for(auto e: solution_edges_in_matching_graph) {
        auto v = matching_to_instance_graph[boost::target(e, matching_graph)];
        auto w = matching_to_instance_graph[boost::source(e, matching_graph)];
        solution.push_back({vertex_to_configuration[v].first, vertex_to_configuration[w].first});
        cost += boost::get(boost::get(boost::edge_weight, matching_graph), e);
    }
    return {cost, solution};
}
