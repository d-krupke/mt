//
// Created by doms on 7/15/16.
//

#include "formulation_3.h"

subset::form3::Solver::Solver(grid_graph& g, std::set<coord>& subset, double turn_cost, double dist_cost, std::time_t deadline)
	: BaseSolver{g, subset, turn_cost, dist_cost, deadline}
{
	std::cout << "Creating IP using formulation 3..." <<std::endl;
	//Build Auxiliary Graph
	for(auto v_it = boost::vertices(g.grid_graph_boost); v_it.first!=v_it.second; ++v_it.first){
		auto v = *v_it.first;
		auto v_north = boost::add_vertex(internal_graph);
		auto v_east = boost::add_vertex(internal_graph);
		auto v_south = boost::add_vertex(internal_graph);
		auto v_west = boost::add_vertex(internal_graph);

		north_map[v] = v_north;
		east_map[v] = v_east;
		south_map[v] = v_south;
		west_map[v] = v_west;
		internal_to_external_graph[v_north] = v;
		internal_to_external_graph[v_east] = v;
		internal_to_external_graph[v_south] = v;
		internal_to_external_graph[v_west] = v;

		boost::add_edge(v_north, v_east, m_TURN_COST, internal_graph);
		boost::add_edge(v_east, v_north, m_TURN_COST, internal_graph);
		boost::add_edge(v_east, v_south, m_TURN_COST, internal_graph);
		boost::add_edge(v_south, v_east, m_TURN_COST, internal_graph);
		boost::add_edge(v_south, v_west, m_TURN_COST, internal_graph);
		boost::add_edge(v_west, v_south, m_TURN_COST, internal_graph);
		boost::add_edge(v_west, v_north, m_TURN_COST, internal_graph);
		boost::add_edge(v_north, v_west, m_TURN_COST, internal_graph);
	}
	for(auto v_it = boost::vertices(g.grid_graph_boost); v_it.first!=v_it.second; ++v_it.first) {
		auto v = *v_it.first;
		for(auto n_it = boost::out_edges(v, g.grid_graph_boost); n_it.first!=n_it.second; ++n_it.first){
			auto e = *n_it.first;
			auto v_1 = boost::target(e, g.grid_graph_boost);
			auto v_2 = boost::source(e, g.grid_graph_boost);
			auto w = (v_1==v?v_2:v_1);

			if(g.getCoord(v).first < g.getCoord(w).first){ //v is west of w
				boost::add_edge(east_map[v], east_map[w], m_DIST_COST, internal_graph);
			} else if(g.getCoord(v).first > g.getCoord(w).first){ // v is east of w
				boost::add_edge(west_map[v], west_map[w], m_DIST_COST, internal_graph);
			} else if(g.getCoord(v).second < g.getCoord(w).second){ //v is south of w
				boost::add_edge(north_map[v], north_map[w], m_DIST_COST, internal_graph);
			} else {
				boost::add_edge(south_map[v], south_map[w], m_DIST_COST, internal_graph);
			}
		}
	}

	//CPLEX

	//Objective
	auto edgeWeightMap = get(boost::edge_weight, internal_graph);
	IloExpr obj_expr{m_env};
	for(auto e_it = boost::edges(internal_graph); e_it.first!=e_it.second; ++e_it.first){
		auto e = *e_it.first;
		auto variable = IloIntVar(m_env, 0, 4);
		variables[e] = variable;
		auto weight = get(edgeWeightMap, e);
		if(weight>0){
			obj_expr+= weight*variable;
		}
	}
	m_objective = IloMinimize(m_env, obj_expr);
	m_model.add(m_objective);

	//enter=leave
	for(auto v_it=boost::vertices(internal_graph); v_it.first!=v_it.second; ++v_it.first){
		auto v = *v_it.first;
		IloExpr constr_expr{m_env};
		for(auto in_it= boost::in_edges(v, internal_graph); in_it.first!=in_it.second; ++in_it.first){
			auto e = *in_it.first;
			constr_expr += 1*variables[e];
		}
		for(auto in_it= boost::out_edges(v, internal_graph); in_it.first!=in_it.second; ++in_it.first){
			auto e = *in_it.first;
			constr_expr += (-1)*variables[e];
		}
		m_model.add(IloRange(m_env, 0, constr_expr, 0));
	}

	//leave>=1 if in subset
	for(auto v_it=boost::vertices(g.grid_graph_boost); v_it.first!=v_it.second; ++v_it.first){
		auto v = *v_it.first;
		if(m_subset.count(m_instance.getCoord(v))>0){
			IloExpr constr_expr{m_env};
			for(auto map: {&north_map, &east_map, &south_map, &west_map}) {
				for (auto e_it = boost::out_edges((*map)[v], internal_graph); e_it.first != e_it.second; ++e_it.first) {
					auto e = *e_it.first;
					if (internal_to_external_graph[boost::source(e, internal_graph)] != internal_to_external_graph[boost::target(e, internal_graph)]) {//external edge
						constr_expr += variables[e];
					}
				}
			}
			m_model.add(IloRange(m_env, 1, constr_expr, 4));
		}
	}

	m_cplex = IloCplex(m_model);
	//if(m_DIST_COST==0) m_cplex.setParam(IloCplex::ObjDif, 0.95*m_TURN_COST);
	std::cout << "Done." <<std::endl;
}


size_t subset::form3::Solver::cut_1(std::map<InternalGraph::edge_descriptor, IloNum>& edge_variable_assignments)
{
	std::cout << "Searching for applications of separation method 1...\n";
	size_t cuts_added = 0;

	//auxiliary graph
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property> comp_graph_t;
	std::map<grid_graph::vertex_descriptor, comp_graph_t::vertex_descriptor> gg_to_cg;
	comp_graph_t cg;

	//add a vertex for each field
	for(auto v_gg_it=boost::vertices(m_instance.grid_graph_boost); v_gg_it.first!=v_gg_it.second; ++v_gg_it.first){
		auto v_gg = *v_gg_it.first;
		auto v_aux = boost::add_vertex(cg);
		gg_to_cg[v_gg] = v_aux;
	}

	//connect all vertices that share a variable>0
	for(auto entry: edge_variable_assignments){
		if(entry.second>0.9) {
			auto uv_aux = entry.first;
			auto u = boost::source(uv_aux, internal_graph);
			auto v = boost::target(uv_aux, internal_graph);
			auto u_gg = internal_to_external_graph[u];
			auto v_gg = internal_to_external_graph[v];
			if(u_gg!=v_gg){
				boost::add_edge(gg_to_cg[u_gg], gg_to_cg[v_gg], cg);
			}
		}
	}

	//Calculate connected components in cg
	std::vector<int> components_map(boost::num_vertices(cg));
	const int num_comps = boost::connected_components(cg, &components_map[0]);
	if(num_comps>1){
		std::map<int, std::set<gg_vd>> components;
		for(auto v_it=boost::vertices(m_instance.grid_graph_boost); v_it.first!=v_it.second; ++v_it.first){
			auto v = *(v_it.first);
			auto v_comp_id = components_map[boost::get(boost::get(boost::vertex_index, cg), gg_to_cg[v])];
			components[v_comp_id].insert(v);
		}

		auto lambda_get_out_edges = [&](gg_vd v)->std::set<InternalGraph::edge_descriptor>{
			std::set<InternalGraph::edge_descriptor> out_edges;
			for(auto v_int: {north_map[v], east_map[v], south_map[v], west_map[v]}){
				for(auto e_p=boost::out_edges(v_int, internal_graph); e_p.first!=e_p.second; ++e_p.first){
					auto w_int = boost::target(*(e_p.first), internal_graph);
					if(internal_to_external_graph[v_int]!=internal_to_external_graph[w_int]){
						out_edges.insert(*(e_p.first));
					}
				}
			}
			return out_edges;
		};

		//create constraints
		for(const auto comp: components){
			size_t card_comp_cap_subset = 0;
			for(auto v: comp.second){
				if(m_subset.count(m_instance.getCoord(v))>0) card_comp_cap_subset+=1;
			}
			if(card_comp_cap_subset>0 && card_comp_cap_subset<m_subset.size()){
				IloExpr constr_expr(m_env);

				for(const auto v: comp.second){
					for(auto e: lambda_get_out_edges(v)){
						constr_expr+= 1*variables[e];
					}
				}

				m_model.add(IloRange(m_env, 1, constr_expr, IloInfinity));
				cuts_added+=1;
			}
		}
	}
	return cuts_added;
}

size_t subset::form3::Solver::eliminate_subtours(subset::SeparationMethod sm)
{
	//obtain variable assignment
	std::map<InternalGraph::edge_descriptor, IloNum> edge_variable_assignments;
	for(auto var_p: variables){
		auto value = m_cplex.getValue(var_p.second);
		edge_variable_assignments[var_p.first]=value;
	}


	if(sm == subset::SeparationMethod::BASIC) return cut_2(edge_variable_assignments);
	return cut_1(edge_variable_assignments)+cut_2(edge_variable_assignments);
}

size_t subset::form3::Solver::cut_2(std::map<InternalGraph::edge_descriptor, IloNum>& solution)
{
	size_t cuts_added = 0;

	//auxilary graph for calculating components (subtours)
	typedef boost::property<boost::vertex_name_t, grid_graph::vertex_descriptor> MapBackToInternalProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, MapBackToInternalProperty, boost::no_property> AuxiliaryGraph;
	AuxiliaryGraph aux_graph;
	std::map<InternalGraph::vertex_descriptor, AuxiliaryGraph::vertex_descriptor> int_to_cg;

	//add vertices to auxilary graph
	for(auto v_it = boost::vertices(internal_graph); v_it.first!=v_it.second;++v_it.first){

		auto v_int = *v_it.first;
		bool v_is_used = false;
		for(auto e_it=boost::out_edges(v_int, internal_graph); e_it.first!=e_it.second; ++e_it.first){
			auto e = *e_it.first;
			if(solution[e]>0.1){
				v_is_used=true;
				break;
			}
		}
		if(v_is_used) {
			auto v_cg = boost::add_vertex(aux_graph);
			auto v_gg = internal_to_external_graph[v_int];
			boost::put(boost::get(boost::vertex_name, aux_graph), v_cg, v_gg);
			int_to_cg[v_int] = v_cg;
		}
	}

	//add edges to auxilary graph
	for(auto e_it = boost::edges(internal_graph); e_it.first!=e_it.second; ++e_it.first){
		auto e_int = *e_it.first;
		if(solution[e_int]>0.1) {
			auto v_int = boost::source(e_int, internal_graph);
			auto w_int = boost::target(e_int, internal_graph);
			boost::add_edge(int_to_cg[v_int], int_to_cg[w_int], aux_graph);
		}
	}

	//calculate components
	std::vector<int> component(boost::num_vertices(aux_graph));
	const auto comp_count = boost::connected_components(aux_graph, &component[0]);
	std::map<int, std::vector<grid_graph::vertex_descriptor>> comps;
	std::map<grid_graph::vertex_descriptor, std::vector<int>> in_comps;
	if(comp_count>1){
		for(auto v_it=boost::vertices(aux_graph); v_it.first!=v_it.second;++v_it.first){
			auto v_cg = *v_it.first;
			auto v_comp_idx = component[boost::get(boost::get(boost::vertex_index, aux_graph), v_cg)];
			auto v_gg = boost::get(boost::get(boost::vertex_name, aux_graph), v_cg);
			if(!comps.count(v_comp_idx)) comps[v_comp_idx]=std::vector<grid_graph::vertex_descriptor>();
			if(std::find(comps[v_comp_idx].begin(), comps[v_comp_idx].end(), v_gg)==comps[v_comp_idx].end()) {
				comps[v_comp_idx].push_back(v_gg);
			}

			if(!in_comps.count(v_gg)) in_comps[v_gg] = std::vector<int>();
			if(std::find(in_comps[v_gg].begin(), in_comps[v_gg].end(), v_comp_idx)==in_comps[v_gg].end()){
				in_comps[v_gg].push_back(v_comp_idx);
			}
		}

		for(auto comp_it: comps){
			assert(comp_it.second.size()>=1);//Actually this should be >1 but './solver.opt -d 900 -f 3 -i ./instances/s_100_sparse_1.gg -t basic -l results_04.csv' failed. 
			//find good field
			gg_vd selected_field;
			size_t rating_p = std::numeric_limits<size_t>::max();
			bool p_is_initialized = false;
			for(auto v: comp_it.second){
				if(m_subset.count(m_instance.getCoord(v))>0){
					size_t rating_v = in_comps[v].size();	
					if(!p_is_initialized || rating_v< rating_p) {
						selected_field = v;
						rating_p = rating_v;
						p_is_initialized = true;
					}
				}
			}

			if(p_is_initialized){
				std::set<gg_vd> cycle;
				for(auto c: in_comps[selected_field]){
					for(auto v: comps[c]){
						cycle.insert(v);
					}
				}

				size_t subset_fields_in_cycle = 0;
				for(auto v: cycle){
					if(m_subset.count(m_instance.getCoord(v))>0) ++subset_fields_in_cycle;
				}
                assert(subset_fields_in_cycle>0);
				if(subset_fields_in_cycle<m_subset.size()) {
					IloExpr constr_expr{m_env};
					for(auto v_ext: cycle) {
						for(auto map: {&north_map, &east_map, &south_map, &west_map}) {
							auto v_int = (*map)[v_ext];
							for(auto in_it=boost::in_edges(v_int, internal_graph); in_it.first!=in_it.second;++in_it.first){
								auto e = *in_it.first;
								if(solution[e]<0.1){
									constr_expr+=1*variables[e];
								}
							}
							for(auto out_it=boost::out_edges(v_int, internal_graph); out_it.first!=out_it.second;++out_it.first){
								auto e = *out_it.first;
								if(solution[e]<0.1){
									constr_expr+=1*variables[e];
								}
							}
						}
					}
					m_model.add(IloRange(m_env, 1, constr_expr, IloInfinity));
					++cuts_added;
				}
			}
		}
		std::cout << "Solution consists of "<<comp_count<<" components. Added "<<cuts_added<< " cuts."<<std::endl;
	} else {
		std::cout << "Solution is tour. No cuts added." <<std::endl;
	}
	assert(comp_count==1 || cuts_added>0);
	return cuts_added;
}
