
#include "formulation_3.h"

penalty::form3::Solver::Solver(grid_graph& g, std::map<coord, double>& penalties, double turn_cost, double dist_cost, std::time_t deadline)
	:BaseSolver{g, penalties, turn_cost, dist_cost, deadline}
{
	std::cout << "Creating IP using formulation 3..." <<std::endl;
	for(auto p: penalties){
		this->penalties_vd[m_instance.getVertex(p.first)] = p.second;
	}

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

		penalty_paid_variables[v] = IloBoolVar{m_env};
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
	for(auto vp: penalty_paid_variables){
		obj_expr+= this->penalties_vd[vp.first]*vp.second;
	}
	IloObjective objective = IloMinimize(m_env, obj_expr);
	m_model.add(objective);

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

	//leave>=1 or penalty paid
	for(auto v_it=boost::vertices(g.grid_graph_boost); v_it.first!=v_it.second; ++v_it.first){
		IloExpr constr_expr{m_env};
		auto v = *v_it.first;
		for(auto map: {&north_map, &east_map, &south_map, &west_map}) {
			for (auto e_it = boost::out_edges((*map)[v], internal_graph); e_it.first != e_it.second; ++e_it.first) {
				auto e = *e_it.first;
				if (internal_to_external_graph[boost::source(e, internal_graph)] != internal_to_external_graph[boost::target(e, internal_graph)]) {//external edge
					constr_expr += variables[e];
				}
			}
		}
		constr_expr += 4*penalty_paid_variables[v];
		m_model.add(IloRange(m_env, 1, constr_expr, 4));
	}

	m_cplex = IloCplex(m_model);
	std::cout << "Done." <<std::endl;
}


size_t penalty::form3::Solver::cut_1(std::map<InternalGraph::edge_descriptor, IloNum>& edge_variable_assignments, std::map<gg_vd, IloNum>& penalty_variable_assigments)
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
			//collect those fields interior and exterior for which the penalty has not been paid
			std::set<gg_vd> active_fields_interior;
			std::set<gg_vd> active_fields_exterior;
			auto lambda_penalty_paid = [&](gg_vd v) -> bool {
				assert(penalty_paid_variables.count(v)>0);
				bool is_true = std::round(static_cast<double>(penalty_variable_assigments[v]))>0;
				return is_true;
			};
			for(auto v: comp.second){
				if(!lambda_penalty_paid(v)) active_fields_interior.insert(v);
			}
			if(active_fields_interior.empty()) continue;
			for(auto v_p = boost::vertices(m_instance.grid_graph_boost); v_p.first!=v_p.second; ++v_p.first){
				auto v = *(v_p.first);
				if(comp.second.count(v)==0 && !lambda_penalty_paid(v)){
					active_fields_exterior.insert(v);
				}
			}
			if(active_fields_exterior.empty()) continue;


			IloExpr constr_expr(m_env);

			for(const auto v: comp.second){
				for(auto e: lambda_get_out_edges(v)){
					constr_expr+= 1*variables[e];
				}
			}

			//penalty part of constraint
			constexpr double REPF = 1.001;//rounding error prevention factor
			for(auto v: active_fields_interior){
				constr_expr+= REPF*(1.0/active_fields_interior.size())*penalty_paid_variables.at(v);
			}
			for(auto v: active_fields_exterior){
				constr_expr+= REPF*(1.0/active_fields_exterior.size())*penalty_paid_variables.at(v);
			}

			m_model.add(IloRange(m_env, 1, constr_expr, IloInfinity));
			cuts_added+=1;
		}
	}
	return cuts_added;
}

size_t penalty::form3::Solver::eliminate_subtours(SeparationMethod sm) {	//solution
	std::map<InternalGraph::edge_descriptor, IloNum> edge_variable_assignments;
	for(auto var_p: variables){
		auto value = m_cplex.getValue(var_p.second.asVariable());
		edge_variable_assignments[var_p.first]=value;
	}
	std::map<gg_vd, IloNum> penalty_variable_assigments;
	for(auto var_p: penalty_paid_variables){
		auto value = m_cplex.getValue(var_p.second);
		penalty_variable_assigments[var_p.first]=value;
	}
	if(sm==SeparationMethod::BASIC){
		return cut_2(edge_variable_assignments, penalty_variable_assigments);
	} else if(sm==SeparationMethod::ADVANCED) {
		return cut_1(edge_variable_assignments, penalty_variable_assigments)+cut_2(edge_variable_assignments, penalty_variable_assigments);
	} else {
        std::cerr << "Illegal eliminate_subtour parameter"<<std::endl;
        std::exit(1);
    }
}

size_t penalty::form3::Solver::cut_2(std::map<InternalGraph::edge_descriptor, IloNum>& edge_variable_assignments, std::map<gg_vd, IloNum>& penalty_variable_assigments)
{
	size_t cuts_added = 0;


	//auxilary graph for calculating components (subtours)
	typedef boost::property<boost::vertex_name_t, grid_graph::vertex_descriptor> MapBackAuxVertexProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, MapBackAuxVertexProperty, boost::no_property> AuxiliaryGraph;
	AuxiliaryGraph aux_graph;
	std::map<InternalGraph::vertex_descriptor, AuxiliaryGraph::vertex_descriptor> int_to_aux;

	//add vertices to auxilary graph
	for(auto v_it = boost::vertices(internal_graph); v_it.first!=v_it.second;++v_it.first){
		auto v_int = *v_it.first;

		//check if vertex is used, otherwise we ignore it
		bool v_is_used = false;
		for(auto e_it=boost::out_edges(v_int, internal_graph); e_it.first!=e_it.second; ++e_it.first){
			auto e = *e_it.first;
			if(edge_variable_assignments[e]>0.1){
				v_is_used=true;
				break;
			}
		}
		//Add vertex to component graph
		if(v_is_used) {
			auto v_cg = boost::add_vertex(aux_graph);
			auto v_gg = internal_to_external_graph[v_int];
			boost::put(boost::get(boost::vertex_name, aux_graph), v_cg, v_gg);
			int_to_aux[v_int] = v_cg;
		}
	}

	//add edges to auxilary graph
	for(auto e_it = boost::edges(internal_graph); e_it.first!=e_it.second; ++e_it.first){
		auto e_int = *e_it.first;
		if(edge_variable_assignments[e_int]>0.1) {
			auto v_int = boost::source(e_int, internal_graph);
			auto w_int = boost::target(e_int, internal_graph);
			auto v_ext = internal_to_external_graph[v_int];// boost::get(boost::get(boost::vertex_name, InternalGraph), v_int).first;
			auto w_ext = internal_to_external_graph[w_int];//boost::get(boost::get(boost::vertex_name, InternalGraph), w_int).first;

			boost::add_edge(int_to_aux[v_int], int_to_aux[w_int], aux_graph);
		}
	}

	//calculate components
	std::vector<int> component(boost::num_vertices(aux_graph));
	int num = boost::connected_components(aux_graph, &component[0]);
	std::map<int, std::set<gg_vd>> comps;
	std::map<gg_vd, std::set<int>> in_comps;
	if(num>1){
		//collect components
		for(auto v_it=boost::vertices(aux_graph); v_it.first!=v_it.second;++v_it.first){
			auto v_cg = *v_it.first;
			auto v_comp_idx = component[boost::get(boost::get(boost::vertex_index, aux_graph), v_cg)];
			auto v_gg = boost::get(boost::get(boost::vertex_name, aux_graph), v_cg);

			comps[v_comp_idx].insert(v_gg);
			in_comps[v_gg].insert(v_comp_idx);
		}

		for(auto comp_it: comps){
			assert(comp_it.second.size()>1);

			//**collect P**
			auto lambda_penalty_paid = [&](gg_vd v) -> bool {
				bool is_true = std::round(static_cast<double>(penalty_variable_assigments[v]))>0;
				return is_true;
			};
			std::vector<gg_vd> P;
			double penalty_sum = 0.0;
			for(auto p: comp_it.second){
				if(in_comps[p].size()==1 && !lambda_penalty_paid(p)) P.push_back(p);
			}
			if(P.empty()){
				gg_vd p = *(comp_it.second.begin());
				auto min_cycle_count = in_comps[p].size();
				for(auto v: comp_it.second){
					auto cycle_count = in_comps[v].size();
					if(cycle_count < min_cycle_count){
						p = v;
						min_cycle_count = cycle_count;
					}
					if(min_cycle_count==1) break;
				}
				P.push_back(p);
			} else {
				//reduce size of P to the expensive fields
				std::sort(P.begin(), P.end(), [&](gg_vd a, gg_vd b)->bool{
						return penalties_vd[a]>penalties_vd[b];
						});
				std::vector<gg_vd> PP;
				double sum=0;
				for(auto p: P){
					if(sum>0.5*penalty_sum && !PP.empty()){
						PP.push_back(p);
						sum+=penalties_vd[p];
					}
				}
				P = PP;
			}

			//collect cycle
			std::set<gg_vd> cycle;

			{
				auto p= *(P.begin());//only more than one comp if |P|==1
				for (auto c: in_comps[p]){
					for(auto v: comps[c]){
						cycle.insert(v);
					}
				}
			}

			std::set<gg_vd> active_fields_exterior;
			for(auto v_p = boost::vertices(m_instance.grid_graph_boost); v_p.first!=v_p.second; ++v_p.first){
				auto v = *(v_p.first);
				if(cycle.count(v)==0 && !lambda_penalty_paid(v)){
					active_fields_exterior.insert(v);
				}
			}
			assert(!active_fields_exterior.empty());

			//add cut: There has to be a change on some field crossed by the subtour (an additional edge has to
			// be used. Removing edges does not help.)
			IloExpr constr_expr{m_env};
			for(auto v_ext: comp_it.second) {
				for(auto map: {&north_map, &east_map, &south_map, &west_map}){
					auto v_int = (*map)[v_ext];
					for(auto in_it=boost::in_edges(v_int, internal_graph); in_it.first!=in_it.second;++in_it.first){
						auto e = *in_it.first;
						if(edge_variable_assignments[e]<0.1){
							constr_expr+=1*variables[e];
						}
					}
					for(auto out_it=boost::out_edges(v_int, internal_graph); out_it.first!=out_it.second;++out_it.first){
						auto e = *out_it.first;
						if(edge_variable_assignments[e]<0.1){
							constr_expr+=1*variables[e];
						}
					}
				}
			}

			constexpr double REPF = 1.001;//rounding error prevention factor
			for(auto v: P){
				constr_expr+= REPF*(1.0/P.size())*penalty_paid_variables[v];
			}
			for(auto v: active_fields_exterior){
				constr_expr+= REPF*(1.0/active_fields_exterior.size())*penalty_paid_variables[v];
			}

			m_model.add(IloRange(m_env, 1, constr_expr, IloInfinity));
			++cuts_added;
		}
        assert(cuts_added>0);
		std::cout << "Solution consists of "<<num<<" components. Added "<<cuts_added<< " cuts."<<std::endl;
	} else {
		std::cout << "Solution is tour. No cuts added." <<std::endl;
	}
	return cuts_added;
}
