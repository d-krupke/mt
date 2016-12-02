
#include "formulation_1.h"
#include <stdexcept>

fc::form1::Solver::Solver(grid_graph& instance, double turn_cost, double dist_cost, std::time_t deadline)
	:BaseSolver{instance, turn_cost, dist_cost, deadline}
{
	std::cout << "Creating IP Formulation for full coverage using the first formulation..." <<std::endl;
	IloExpr obj_expr(m_env);
	std::map<grid_graph::boost_grid_graph_t::edge_descriptor, IloExpr> edge_constr_map; // every edge has to be passed the same amount from each side
	std::map<grid_graph::boost_grid_graph_t::vertex_descriptor, IloExpr> vertex_constr_map; //every vertex has to be visited at least once

	for(auto v_it_p=boost::vertices(m_instance.grid_graph_boost); v_it_p.first!=v_it_p.second; ++v_it_p.first){
		const auto v = *v_it_p.first;
		for(auto w_it_p=boost::adjacent_vertices(v, m_instance.grid_graph_boost); w_it_p.first!=w_it_p.second; ++w_it_p.first) {
			const auto w = *w_it_p.first; assert(v!=w);
            assert(v != w);

			//u-turn w-v-w
			auto var_wvw = IloIntVar(m_env,0, 4, (std::to_string(v)+"_"+std::to_string(w)).c_str());
			variable_map[var_index{w,v,w}] = var_wvw;
			obj_expr+=(2*m_TURN_COST+m_DIST_COST)*var_wvw;
            //wvw==True -> v is visited
			if(!vertex_constr_map.count(v)) vertex_constr_map[v] = IloExpr(m_env);
			vertex_constr_map[v] += 1*var_wvw;
			const auto vw = boost::edge(v,w, m_instance.grid_graph_boost).first;
			if(!edge_constr_map.count(vw)) edge_constr_map[vw] = IloExpr(m_env);
			edge_constr_map[vw]+= (v<w?2:-2)*var_wvw;

			//other turns u-v-w with u!=w
			for (auto u_it_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); u_it_p.first != u_it_p.second; ++u_it_p.first) {
				const auto u = *u_it_p.first;
				if(u<=w) continue;

				auto var_uvw = IloIntVar(m_env,0, 4,(std::to_string(u)+"_"+std::to_string(v)+"_"+std::to_string(w)).c_str());
				assert(!variable_map.count(var_index{u,v,w}));
				variable_map[var_index{u,v,w}] = var_uvw;

                //uvw==True -> v is visited
				if(!vertex_constr_map.count(v)) vertex_constr_map[v] = IloExpr(m_env);
				vertex_constr_map[v] += 1*var_uvw;

				auto uv = boost::edge(u,v, m_instance.grid_graph_boost).first;
				auto vw = boost::edge(v,w, m_instance.grid_graph_boost).first;
				if(!edge_constr_map.count(uv)) edge_constr_map[uv] = IloExpr(m_env);
				if(!edge_constr_map.count(vw)) edge_constr_map[vw] = IloExpr(m_env);
				edge_constr_map[uv]+= (v<u?1:-1)*var_uvw;
				edge_constr_map[vw]+= (v<w?1:-1)*var_uvw;

				const auto u_pos = m_instance.getCoord(u);
				const auto w_pos = m_instance.getCoord(w);
				if((w_pos.first != u_pos.first) && (w_pos.second != u_pos.second)){ //90deg turn
					obj_expr+=(m_TURN_COST+m_DIST_COST)*var_uvw;
				} else { //straight pass
					obj_expr+=m_DIST_COST*var_uvw;
				}
			}
		}
	}

	//Add objective and constraints to model
	m_objective = IloMinimize(m_env, obj_expr);
	m_model.add(m_objective);

	for(auto e_expr_it: edge_constr_map){
		IloRange constraint(m_env, 0, e_expr_it.second, 0); //in - out = 0
		m_model.add(constraint);
	}

	for(auto v_expr_it: vertex_constr_map){
		IloRange constraint(m_env, 1, v_expr_it.second, 4); //At least 1 visit per field
		m_model.add(constraint);
	}

	//init Solver
	m_cplex = IloCplex(m_model);
	//if(m_DIST_COST==0) m_cplex.setParam(IloCplex::ObjDif, 0.95*m_TURN_COST);

	std::cout << "Created IP Formulation." <<std::endl;
}

size_t fc::form1::Solver::eliminate_subtours(fc::SeparationMethod sm)
{
	if(sm == fc::SeparationMethod::NONE) throw std::invalid_argument{"Eliminating subtours with NONE"};
	std::map<var_index, IloNum> edge_variable_assignments;
	for(auto p: variable_map){
		IloNum v = m_cplex.getValue(p.second);
		edge_variable_assignments[p.first] = v;
	}
	if(sm == fc::SeparationMethod::BASIC) return cut_2(edge_variable_assignments);
	return cut_1(edge_variable_assignments)+cut_2(edge_variable_assignments);
}

size_t fc::form1::Solver::cut_1(std::map<var_index, IloNum>& edge_variable_assignments)
{
	std::cout << "Searching for applications of simple insufficient separation constraint...\n";
	size_t cuts_added = 0;

	//auxiliary graph
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property> aux_graph_t;
	std::map<grid_graph::vertex_descriptor, aux_graph_t::vertex_descriptor> gg_to_aux;
	aux_graph_t aux_graph;

	//add a vertex for each field
	for(auto v_gg_it=boost::vertices(m_instance.grid_graph_boost); v_gg_it.first!=v_gg_it.second; ++v_gg_it.first){
		auto v_gg = *v_gg_it.first;
		auto v_aux = boost::add_vertex(aux_graph);
		gg_to_aux[v_gg] = v_aux;
	}

	//connect all vertices that share a variable>0
	for(auto entry: edge_variable_assignments){
		if(entry.second>0.9) {
			auto u_gg = std::get<0>(entry.first);
			auto v_gg = std::get<1>(entry.first);
			auto w_gg = std::get<2>(entry.first);

			boost::add_edge(gg_to_aux[u_gg], gg_to_aux[v_gg], aux_graph);
			boost::add_edge(gg_to_aux[v_gg], gg_to_aux[w_gg], aux_graph);
		}
	}

	//Calculate connected components in InternalGraph
	std::vector<int> components_map(boost::num_vertices(aux_graph));
	const int num_comps = boost::connected_components(aux_graph, &components_map[0]);
	if(num_comps>1){
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

		//create constraints
		for(const auto comp: components){
			IloExpr constr_expr(m_env);

			for(const auto v: comp.second){
				for(const auto uvw: lambda_get_variable_indexes(v)){
					const auto u = std::get<0>(uvw);
					const auto w = std::get<2>(uvw);
					if(comp.second.count(u)==0 || comp.second.count(w)==0){
						constr_expr+=1*variable_map[uvw];
					}
				}
			}

			m_model.add(IloRange(m_env, 1, constr_expr, IloInfinity));
			cuts_added+=1;

            std::cout << "Added "<< cuts_added << " new constraints.\n";
		}
	} else {
        std::cout << "No application of this separation possible."<< std::endl;
    }
	return cuts_added;
}

size_t fc::form1::Solver::cut_2(std::map<var_index, IloNum>& edge_variable_assignments)
{
	std::cout << "Searching for applications of advanced sufficient separation constraint..." << std::endl;
	size_t cuts_added = 0;

    /**
     * Calculating the components
     * ----------------------------------------------
     * We use an auxiliary graph in which every edge in the grid graph is a vertex.
     * Two vertices e.g. (u,v) and (v,w) are connected if there is a passing that ends in both, i.e., the passing uvw (covering v and ending at u and w).
     * If and only if two cycles share a same vertex in this auxiliary graph, they can be connected for free.
     * Calculating the maximal cycles hence becomes calculating the connected components that contain an edge in the auxiliary graph.
     * Of course there will be lots of isolated vertices for which we do not care.
     */
	struct cycle_cut_graph_idx: std::pair<gg_vd,gg_vd> //Index for the vertices in the auxiliary graph
    {
		cycle_cut_graph_idx(): std::pair<gg_vd,gg_vd>{}{ /*-*/ }
		cycle_cut_graph_idx(gg_vd a, gg_vd b): std::pair<gg_vd,gg_vd>{std::max(a,b), std::min(a,b)}{ /*-*/ }
	};

    //Defining the auxiliary graph
	typedef boost::property<boost::vertex_name_t, cycle_cut_graph_idx> vertex_name_property_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, vertex_name_property_t> boost_cycle_cut_graph_t;
	boost_cycle_cut_graph_t aux_graph;
	std::map<cycle_cut_graph_idx, boost_cycle_cut_graph_t::vertex_descriptor > vertex_map; //Get the boost vertex from {u,w}
	std::map<gg_vd, std::vector<boost_cycle_cut_graph_t::vertex_descriptor>> field_to_aux_vertices; //Get all vertices for a field (up to four)

    //Building the auxiliary graph
	for(auto va: edge_variable_assignments){
		if(va.second>0.9) { // loop over all used passings uvw
			auto u = std::get<0>(va.first);
			auto v = std::get<1>(va.first);
			auto w = std::get<2>(va.first);
			cycle_cut_graph_idx uv{u,v};
			cycle_cut_graph_idx vw{v,w};
			if(!vertex_map.count(uv)){
				vertex_map[uv] = boost::add_vertex(aux_graph);
				field_to_aux_vertices[u].push_back(vertex_map[uv]);
				field_to_aux_vertices[v].push_back(vertex_map[uv]);
			}
			if(!vertex_map.count(vw)){
				vertex_map[vw] = boost::add_vertex(aux_graph);
				field_to_aux_vertices[v].push_back(vertex_map[vw]);
				field_to_aux_vertices[w].push_back(vertex_map[vw]);
			}
			if(uv!=vw){
				boost::add_edge(vertex_map[uv], vertex_map[vw], aux_graph);
			}
		}
	}


	//Calculate connected components
	std::vector<int> component(boost::num_vertices(aux_graph));
	const auto num_comps = boost::connected_components(aux_graph, &component[0]);

	if(num_comps>1){//if there is more than one connected component we have to separate

		// ----------------- auxiliary functions ---------------------------
		auto is_straight_field = [&](gg_vd v){
			for(auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first!=it_u_p.second; ++it_u_p.first){
				for(auto it_w_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_w_p.first!=it_w_p.second; ++it_w_p.first){
					auto u = *it_u_p.first;
					auto w = *it_w_p.first;
					auto pos_u = m_instance.getCoord(u);
					auto pos_w = m_instance.getCoord(w);
					if(pos_u == pos_w || (pos_u.first!=pos_w.first && pos_w.second!=pos_u.second)){
						if(edge_variable_assignments[var_index{u,v,w}]>0.9) return false;
					}
				}
			}
			return true;
		};
		auto get_unused_variables_excluding_u_turns = [&](gg_vd v)->std::set<var_index>{
			std::set<var_index> ret;
			for(auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first!=it_u_p.second; ++it_u_p.first) {
				for (auto it_w_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_w_p.first != it_w_p.second; ++it_w_p.first) {
					auto u = *it_u_p.first;
					auto w = *it_w_p.first;
					if(edge_variable_assignments[var_index{u,v,w}]<0.1){
						ret.insert(var_index{u,v,w});
					}
				}
			}
			return ret;
		};
		auto get_turn_variables = [&](gg_vd v)->std::set<var_index>{
			std::set<var_index> ret;
			for(auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first!=it_u_p.second; ++it_u_p.first) {
				for (auto it_w_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_w_p.first != it_w_p.second; ++it_w_p.first) {
					auto u = *it_u_p.first;
					auto w = *it_w_p.first;
					auto pos_u = m_instance.getCoord(*it_u_p.first);
					auto pos_w = m_instance.getCoord(*it_w_p.first);
					if(pos_u != pos_w && (pos_u.first!=pos_w.first && pos_w.second!=pos_u.second)){
						ret.insert(var_index{u,v,w});
					}
				}
			}
			return ret;
		};
		auto get_components= [&](gg_vd v)->std::set<int>{
			std::set<int> comps;
			auto index_map = boost::get(boost::vertex_index, aux_graph);
			for(auto x: field_to_aux_vertices[v]){
				comps.insert(component[boost::get(index_map, x)]);
			}
			return comps;
		};

		// ------ Build cuts -----------------
		std::map<int, IloExpr> comp_to_expr_map;
		std::map<int, bool> uniquely_visited_field_added;

		for(auto it_v_p = boost::vertices(m_instance.grid_graph_boost); it_v_p.first!=it_v_p.second; ++it_v_p.first){
			auto v = *it_v_p.first;
			auto cycles = get_components(v);
			bool is_unique_visited = cycles.size()==1;
			bool is_straight = is_straight_field(v);
			for(auto c: cycles){
				if(!comp_to_expr_map.count(c)){
					comp_to_expr_map[c] = IloExpr(m_env);
					uniquely_visited_field_added[c] = false;
				}
				if(is_unique_visited && !uniquely_visited_field_added[c]){
					//Add unused variables
					for(auto var: get_unused_variables_excluding_u_turns(v)){ //all unused variables except u-turns
						comp_to_expr_map[c]+=1*variable_map[var];
						assert(edge_variable_assignments[var]<0.1);
					}
					//unused u-turns
					for(auto it_u_p = boost::adjacent_vertices(v, m_instance.grid_graph_boost); it_u_p.first!=it_u_p.second; ++it_u_p.first) {
						auto u = *it_u_p.first;
						if(edge_variable_assignments[var_index{u,v,u}]<0.1) {
							comp_to_expr_map[c]+=variable_map[var_index{u,v,u}];
							auto var_i = var_index{u,v,u};
							assert(edge_variable_assignments[var_i]<0.1);
						}
					}
					uniquely_visited_field_added[c] = true;
				} else if(is_straight){
					for(auto var: get_turn_variables(v)){ //all unused variables except u-turns
						comp_to_expr_map[c]+=1*variable_map[var];
						assert(edge_variable_assignments[var]<0.1);
					}
				} else {
					for(auto var: get_unused_variables_excluding_u_turns(v)){ //all unused variables except u-turns
						comp_to_expr_map[c]+=1*variable_map[var];
						assert(edge_variable_assignments[var]<0.1);
					}
				}
			}
		}

		//--- add cuts ----------------
		for(auto it: comp_to_expr_map){
			if(uniquely_visited_field_added[it.first]){
				IloRange constr{m_env, 1, it.second, IloInfinity};
				m_model.add(constr);
				cuts_added+=1;
			}
		}
		assert(cuts_added<1000); //this should be a rational number
	}
	if(cuts_added>0) std::cout << "Found and added "<<cuts_added<<" separation constraints!"<<std::endl;
	else std::cout << "No cycle cuts found, the cycle cover is a tour!" << std::endl;
	return cuts_added;
}


