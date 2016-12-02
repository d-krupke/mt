
#include "formulation_1.h"

penalty::form1::Solver::Solver(grid_graph &instance, std::map<coord, double>& penalties, double turn_cost, double dist_cost, std::time_t deadline):
	BaseSolver{instance, penalties, turn_cost, dist_cost, deadline}
{
	std::cout << "Creating IP Formulation for penalty coverage using the first formulation..." << std::endl;
	for(auto p: penalties){
		m_penalties_vd[m_instance.getVertex(p.first)] = p.second;
	}

	IloExpr obj_expr(m_env);
	std::map<grid_graph::boost_grid_graph_t::edge_descriptor, IloExpr> edge_constr_map; // every edge has to be passed the same amount from each side

	for (auto v_it_p = boost::vertices(m_instance.grid_graph_boost); v_it_p.first != v_it_p.second; ++v_it_p.first) {
		auto v = *v_it_p.first;
		penalty_paid_variable_map[v]=IloBoolVar(m_env);
		IloExpr visited_expr{m_env};

		for (auto w_it_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); w_it_p.first != w_it_p.second; ++w_it_p.first) {
			auto w = *w_it_p.first;
			assert(v != w);

			// U-Turn (w,v,w)
            {
                auto wvw_var = IloIntVar(m_env, 0, 4, (std::to_string(v) + "_" + std::to_string(w)).c_str());
                variable_map[var_index{w, v, w}] = wvw_var;
                obj_expr += (2 * m_TURN_COST + m_DIST_COST) * wvw_var;

                visited_expr += 1 * wvw_var;

                auto vw = boost::edge(v, w, m_instance.grid_graph_boost).first;
                if (!edge_constr_map.count(vw)) edge_constr_map[vw] = IloExpr(m_env);
                edge_constr_map[vw] += (v < w ? 2 : -2) * wvw_var;
            }

			for (auto u_it_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); u_it_p.first != u_it_p.second; ++u_it_p.first) {
				auto u = *u_it_p.first;
				if (u <= w) continue;

				// Straight or 90deg turn (u,v,w)
				auto uvw_var = IloIntVar(m_env, 0, 4, (std::to_string(u) + "_" + std::to_string(v) + "_" + std::to_string(w)).c_str());
				assert(!variable_map.count(var_index{u, v, w}));
				variable_map[var_index{u, v, w}] = uvw_var;

				visited_expr += 1 * uvw_var;

				auto uv = boost::edge(u, v, m_instance.grid_graph_boost).first;
				auto vw = boost::edge(v, w, m_instance.grid_graph_boost).first;
				if (!edge_constr_map.count(uv)) edge_constr_map[uv] = IloExpr(m_env);
				if (!edge_constr_map.count(vw)) edge_constr_map[vw] = IloExpr(m_env);
				edge_constr_map[uv] += (v < u ? 1 : -1) * uvw_var;
				edge_constr_map[vw] += (v < w ? 1 : -1) * uvw_var;

				const auto u_pos = m_instance.getCoord(u);
				const auto w_pos = m_instance.getCoord(w);
                if (w_pos.first != u_pos.first && w_pos.second != u_pos.second) { //90deg turn
                    obj_expr += (m_TURN_COST + m_DIST_COST) * uvw_var;
                } else { //straight pass
                    obj_expr += (m_DIST_COST) * uvw_var;
                } //u-turns are already added earlier
			}
		}
		//Penalty if not used
		obj_expr += m_penalties_vd[v]*penalty_paid_variable_map.at(v);

		m_model.add((visited_expr+penalty_paid_variable_map.at(v))>=1);
		m_model.add((visited_expr+4*penalty_paid_variable_map.at(v))<=4);
	}
	m_objective = IloMinimize(m_env, obj_expr);
	m_model.add(m_objective);

	for (auto e_expr_it: edge_constr_map) {
		IloRange constraint(m_env, 0, e_expr_it.second, 0); //in - out = 0
		m_model.add(constraint);
	}

	m_cplex = IloCplex(m_model);

	std::cout << "Created IP Formulation." << std::endl;
}


size_t penalty::form1::Solver::eliminate_subtours(SeparationMethod sm)
{
	//variable assignment
	std::map<var_index, IloNum> edge_variable_assignments;
	for (auto p: variable_map) {
		IloNum value = m_cplex.getValue(p.second);
		edge_variable_assignments[p.first] = value;
	}
	std::map<gg_vd, IloNum> penalty_variable_assigments;
	for(auto p: penalty_paid_variable_map){
		IloNum value = m_cplex.getValue(p.second);
		penalty_variable_assigments[p.first]=value;
	}

	if(sm == SeparationMethod ::BASIC){
		return cut_2(edge_variable_assignments, penalty_variable_assigments);
	} else if(sm == SeparationMethod::ADVANCED){
		return cut_1(edge_variable_assignments, penalty_variable_assigments)+cut_2(edge_variable_assignments, penalty_variable_assigments);
	} else {
        std::cerr << "Eliminate subtours with SeparationMethod::NONE"<<std::endl;
        std::exit(1);
    }
}

size_t penalty::form1::Solver::cut_1(std::map<var_index, IloNum>& edge_variable_assignments, std::map<gg_vd, IloNum>& penalty_variable_assigments)
{
	std::cout << "Searching for applications of separation method 1...\n";
	size_t cuts_added = 0;



	//auxiliary graph
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property> aux_graph_t;
	std::map<grid_graph::vertex_descriptor, aux_graph_t::vertex_descriptor> gg_to_aux;
	aux_graph_t aux_graph;

	//add a vertex for each field
	for (auto v_gg_it = boost::vertices(m_instance.grid_graph_boost); v_gg_it.first != v_gg_it.second; ++v_gg_it.first) {
		auto v_gg = *v_gg_it.first;
		auto v_aux = boost::add_vertex(aux_graph);
		gg_to_aux[v_gg] = v_aux;
	}

	//connect all vertices that share a variable>0
	for (auto entry: edge_variable_assignments) {
		if (entry.second > 0.9) {
			auto u_gg = std::get<0>(entry.first);
			auto v_gg = std::get<1>(entry.first);
			auto w_gg = std::get<2>(entry.first);

			boost::add_edge(gg_to_aux[u_gg], gg_to_aux[v_gg], aux_graph);
			boost::add_edge(gg_to_aux[v_gg], gg_to_aux[w_gg], aux_graph);
		}
	}



	//Calculate connected components in InternalGraph
	std::vector<int> components_map(boost::num_vertices(aux_graph));
	int num_comps = boost::connected_components(aux_graph, &components_map[0]);
	std::map<int, std::set<gg_vd>> components;
	for(auto v_it=boost::vertices(m_instance.grid_graph_boost); v_it.first!=v_it.second; ++v_it.first){
		auto v = *(v_it.first);
		auto v_comp_id = components_map[boost::get(boost::get(boost::vertex_index, aux_graph), gg_to_aux[v])];
		components[v_comp_id].insert(v);
	}

	//auxiliary function
	auto lambda_get_variable_indexes=[&](gg_vd v)->std::set<var_index>{
		std::set<var_index> ret;
		for(auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first!=it_u_p.second; ++it_u_p.first) {
			for (auto it_w_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_w_p.first != it_w_p.second; ++it_w_p.first) {
				auto u = *it_u_p.first;
				auto w = *it_w_p.first;
				ret.insert(var_index{u,v,w});
			}
		}
		return ret;
	};

	//go over all components and possibly apply separation. Note that there will be many components of size one without a cycle (penalty paid)
	for(const auto comp: components){
		//collect those fields interior and exterior for which the penalty has not been paid
		std::set<gg_vd> active_fields_interior;
		std::set<gg_vd> active_fields_exterior;
		auto lambda_penalty_paid = [&](gg_vd v) -> bool {
			assert(penalty_paid_variable_map.count(v)>0);
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

		//create constraints
		IloExpr constr_expr(m_env);

		//basic part of constraint
		for(const auto v: comp.second){
			for(const auto uvw: lambda_get_variable_indexes(v)){
				const auto u = std::get<0>(uvw);
				const auto w = std::get<2>(uvw);
				if(comp.second.count(u)==0 || comp.second.count(w)==0){
					constr_expr+=1*variable_map.at(uvw);
				}
			}
		}

		//penalty part of constraint
		constexpr double REPF = 1.001;//rounding error prevention factor
		for(auto v: active_fields_interior){
			constr_expr+= REPF*(1.0/active_fields_interior.size())*penalty_paid_variable_map.at(v);
		}
		for(auto v: active_fields_exterior){
			constr_expr+= REPF*(1.0/active_fields_exterior.size())*penalty_paid_variable_map.at(v);
		}

		m_model.add(IloRange(m_env, 1, constr_expr, IloInfinity));
		cuts_added+=1;
	}

	std::cout << "Added "<< cuts_added << " new constraints.\n";
	return cuts_added;
}

size_t penalty::form1::Solver::cut_2(std::map<var_index, IloNum>& edge_variable_assignments, std::map<gg_vd, IloNum>& penalty_variable_assigments)
{
	std::cout << "Searching for cycle cuts..." << std::endl;
	size_t cuts_added = 0;



	// **** Auxiliary graph *****
	struct cycle_cut_graph_idx : std::pair<gg_vd, gg_vd> {
		cycle_cut_graph_idx() : std::pair<gg_vd, gg_vd>{} {}

		cycle_cut_graph_idx(gg_vd a, gg_vd b) : std::pair<gg_vd, gg_vd>{std::max(a, b), std::min(a, b)} {}
	};

	typedef boost::property<boost::vertex_name_t, cycle_cut_graph_idx> vertex_name_property_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, vertex_name_property_t> boost_cycle_cut_graph_t;

	boost_cycle_cut_graph_t aux_graph;
	std::map<cycle_cut_graph_idx, boost_cycle_cut_graph_t::vertex_descriptor> vertex_map;


	std::map<gg_vd, std::vector<boost_cycle_cut_graph_t::vertex_descriptor>> field_to_aux_vertices;
	for (auto va: edge_variable_assignments) {
		if (va.second > 0.9) {
			auto u = std::get<0>(va.first);
			auto v = std::get<1>(va.first);
			auto w = std::get<2>(va.first);
			if (!field_to_aux_vertices.count(u)) field_to_aux_vertices[u] = std::vector<gg_vd>();
			if (!field_to_aux_vertices.count(v)) field_to_aux_vertices[v] = std::vector<gg_vd>();
			if (!field_to_aux_vertices.count(w)) field_to_aux_vertices[w] = std::vector<gg_vd>();
			cycle_cut_graph_idx uv{u, v};
			cycle_cut_graph_idx vw{v, w};
			if (!vertex_map.count(uv)) {
				vertex_map[uv] = boost::add_vertex(aux_graph);
				field_to_aux_vertices[u].push_back(vertex_map[uv]);
				field_to_aux_vertices[v].push_back(vertex_map[uv]);
			}
			if (!vertex_map.count(vw)) {
				vertex_map[vw] = boost::add_vertex(aux_graph);
				field_to_aux_vertices[v].push_back(vertex_map[vw]);
				field_to_aux_vertices[w].push_back(vertex_map[vw]);
			}
			if (uv != vw) {
				boost::add_edge(vertex_map[uv], vertex_map[vw], aux_graph);
			}
		}
	}


	//Calculate connected components components
	std::vector<int> component(boost::num_vertices(aux_graph));
	int num_comps = boost::connected_components(aux_graph, &component[0]);

	if (num_comps > 1) {//if there is more than one connected component we have to separate

		// ----------------- auxiliary functions ---------------------------
		auto is_straight_field = [&](gg_vd v) {
			for (auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first != it_u_p.second; ++it_u_p.first) {
				for (auto it_w_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_w_p.first != it_w_p.second; ++it_w_p.first) {
					auto u = *it_u_p.first;
					auto w = *it_w_p.first;
					auto pos_u = m_instance.getCoord(u);
					auto pos_w = m_instance.getCoord(w);
					if (pos_u == pos_w || (pos_u.first != pos_w.first && pos_w.second != pos_u.second)) {
						if (edge_variable_assignments[var_index{u, v, w}] > 0.9) return false;
					}
				}
			}
			return true;
		};
		auto get_unused_variables_excluding_u_turns = [&](gg_vd v) -> std::set<var_index> {
			std::set<var_index> ret;
			for (auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first != it_u_p.second; ++it_u_p.first) {
				for (auto it_w_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_w_p.first != it_w_p.second; ++it_w_p.first) {
					auto u = *it_u_p.first;
					auto w = *it_w_p.first;
					if (edge_variable_assignments[var_index{u, v, w}] < 0.1) {
						ret.insert(var_index{u, v, w});
					}
				}
			}
			return ret;
		};
		auto get_turn_variables = [&](gg_vd v) -> std::set<var_index> {
			std::set<var_index> ret;
			for (auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first != it_u_p.second; ++it_u_p.first) {
				for (auto it_w_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_w_p.first != it_w_p.second; ++it_w_p.first) {
					auto u = *it_u_p.first;
					auto w = *it_w_p.first;
					auto pos_u = m_instance.getCoord(*it_u_p.first);
					auto pos_w = m_instance.getCoord(*it_w_p.first);
					if (pos_u != pos_w && (pos_u.first != pos_w.first && pos_w.second != pos_u.second)) {
						ret.insert(var_index{u, v, w});
					}
				}
			}
			return ret;
		};
		auto get_components = [&](gg_vd v) -> std::set<int> {
			std::set<int> comps;
			auto index_map = boost::get(boost::vertex_index, aux_graph);
			for (auto x: field_to_aux_vertices[v]) {
				comps.insert(component[boost::get(index_map, x)]);
			}
			return comps;
		};
		auto lambda_penalty_paid = [&](gg_vd v) -> bool {
			assert(penalty_paid_variable_map.count(v)>0);
			bool is_true = std::round(static_cast<double>(penalty_variable_assigments[v]))>0;
			return is_true;
		};

		std::map<int, std::set<gg_vd>> vertex_components;
		for (int i = 0; i < num_comps; ++i) { vertex_components[i] = std::set<gg_vd>(); }
		for (auto it_v_p = boost::vertices(m_instance.grid_graph_boost); it_v_p.first != it_v_p.second; ++it_v_p.first) {
			auto v = *it_v_p.first;
			auto cycles = get_components(v);
			for (auto c_id: cycles) {
				vertex_components[c_id].insert(v);
			}
		}
		for (auto cycle_entry: vertex_components) {
			assert(cycle_entry.second.size()>1);

			//**Collect a good set P**
			double penalty_sum = 0.0;
			std::vector<gg_vd> P;
			for(auto v: cycle_entry.second){
				auto cycle_count = get_components(v).size();
				if(cycle_count==1 && !lambda_penalty_paid(v)){
					P.push_back(v);
					penalty_sum+=m_penalties_vd[v];
				}
			}
			if(P.empty()){
				//find a good p
				gg_vd p = *cycle_entry.second.begin();
				auto min_cycle_count = get_components(p).size();
				for (auto v: cycle_entry.second) {
					auto cycle_count = get_components(v).size();
					if (cycle_count < min_cycle_count) {
						p = v;
						min_cycle_count = cycle_count;
					}
					if (min_cycle_count == 1) break;
				}
				P.push_back(p);
			} else {
				//reduce size of P to the expensive fields
				std::sort(P.begin(), P.end(), [&](gg_vd a, gg_vd b)->bool{
						return m_penalties_vd[a]>m_penalties_vd[b];
						});
				std::vector<gg_vd> PP;
				double sum=0;
				for(auto p: P){
					if(sum>0.5*penalty_sum && !PP.empty()){
						PP.push_back(p);
						sum+=m_penalties_vd[p];
					}
				}
				P = PP;
			}



			//**collect fields in cycle (possibly consists of multiple subcycles)**
			std::set<gg_vd> cycle;
			for (auto c: get_components(*P.begin())) {//only more than one comp if |P|==1
				cycle.insert(vertex_components[c].begin(), vertex_components[c].end());
			}

			//**create constraints**
			IloExpr constr_expr(m_env);
			for (auto v: cycle) {
				if (std::find(P.begin(), P.end(), v) != P.end()) {//For all fields in P we have to accept any change
					//Add unused variables
					for (auto var: get_unused_variables_excluding_u_turns(v)) { //all unused variables except u-turns
						constr_expr += 1 * variable_map.at(var);
						assert(edge_variable_assignments[var] < 0.1);
					}
					//unused u-turns
					for (auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first != it_u_p.second; ++it_u_p.first) {
						auto u = *it_u_p.first;
						if (edge_variable_assignments[var_index{u, v, u}] < 0.1) {
							constr_expr += variable_map.at(var_index{u, v, u});
							auto var_i = var_index{u, v, u};
							assert(edge_variable_assignments[var_i] < 0.1);
						}
					}

				} else { //for the rest of the fields, we can be more precise.
					if (is_straight_field(v)) {
						for (auto var: get_turn_variables(v)) { //all unused variables except u-turns
							constr_expr += 1 * variable_map.at(var);
							assert(edge_variable_assignments[var] < 0.1);
						}
					} else {
						for (auto var: get_unused_variables_excluding_u_turns(v)) { //all unused variables except u-turns
							constr_expr += 1 * variable_map.at(var);
							assert(edge_variable_assignments[var] < 0.1);
						}
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

			constexpr double REPF = 1.001;//rounding error prevention factor
			for(auto v: P){
				constr_expr+= REPF*(1.0/P.size())*penalty_paid_variable_map.at(v);
			}
			for(auto v: active_fields_exterior){
				constr_expr+= REPF*(1.0/active_fields_exterior.size())*penalty_paid_variable_map.at(v);
			}

			//add constr
			m_model.add(IloRange(m_env, 1, constr_expr, IloInfinity));
			cuts_added += 1;
		}
		assert(cuts_added < 1000); //this should be a rational number
        assert(cuts_added>0);
	}
	if (cuts_added > 0) std::cout << "Found and added " << cuts_added << " cycle cuts!" << std::endl;
	else std::cout << "No cycle cuts found, the cycle cover is a tour!" << std::endl;
	return cuts_added;
}

