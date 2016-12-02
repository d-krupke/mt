//
// Created by doms on 6/24/16.
//

#include "formulation_2.h"
#include <stdexcept>

subset::form2::Solver::Solver(grid_graph& instance, std::set<coord>& subset, double turn_cost, double dist_cost, std::time_t deadline, bool cycle_cover):
	BaseSolver{instance, subset, turn_cost, dist_cost, deadline}	
{
	std::cout << "Creating solver using formulation 2..."<<std::endl;
	IloExpr obj_expr(m_env);
	distance_calculator dc{m_instance, m_TURN_COST, m_DIST_COST}; //For calculating the weights for the internal graph

	// create vertices
	for(auto it= boost::vertices(m_instance.grid_graph_boost); it.first!=it.second; ++it.first){
		auto v = *it.first;
		if(m_subset.count(m_instance.getCoord(v))>0){
			subset_vertices.push_back(v);
			north_map[v] = boost::add_vertex(graph);
			east_map[v] = boost::add_vertex(graph);
			south_map[v] = boost::add_vertex(graph);
			west_map[v] = boost::add_vertex(graph);
			internal_to_grid_graph[north_map[v]] = v;
			internal_to_grid_graph[east_map[v]] = v;
			internal_to_grid_graph[south_map[v]] = v;
			internal_to_grid_graph[west_map[v]] = v;
		}
	}
	std::cout << "subset_vertices.size()="<<subset_vertices.size()<<std::endl;

	//create edges
	for(auto v: subset_vertices){
		for(auto w: subset_vertices){
			if(v>=w) continue;

			auto auto_map = [&](grid_graph::vertex_descriptor vertex, directions d){
				switch(d){
					case directions::north: return north_map[vertex];
					case directions::east: return east_map[vertex];
					case directions::south: return south_map[vertex];
					default: return west_map[vertex];
				}
			};

			for(auto d1: {directions::north,directions::east,directions::south,directions::west}){
				for(auto d2: {directions::north,directions::east,directions::south,directions::west}){
					auto x = boost::add_edge(auto_map(v, d1), auto_map(w, d2), dc.distance(v, d1, w, d2), graph);
					variables[x.first] = IloBoolVar(m_env);
					obj_expr+=(variables[x.first].asVariable()*dc.distance(v, d1, w, d2));
				}
			}
		}
		//Add single field cycles
		if(cycle_cover){
			auto single_field_cycle_cost = 4*turn_cost+2*dist_cost;
			auto e = boost::add_edge(west_map[v], east_map[v], single_field_cycle_cost, graph);
			variables[e.first] = IloBoolVar(m_env);
			obj_expr+=variables[e.first].asVariable()*single_field_cycle_cost;
		}
	}
	IloObjective objective = IloMinimize(m_env, obj_expr);
	m_model.add(objective);

	//add constraint: same amount of edges at each side used
	for(auto v: subset_vertices){
		IloExpr visited_constr{m_env};

		//north-south
		IloExpr nw_constr{m_env};
		for(auto e_p=boost::out_edges(north_map[v], graph); e_p.first!=e_p.second; ++e_p.first){
			auto e = *(e_p.first);
			auto e_var = variables[e];
			nw_constr+=1*e_var;
			visited_constr+=e_var;
		}
		for(auto e_p=boost::out_edges(south_map[v], graph); e_p.first!=e_p.second; ++e_p.first){
			auto e = *(e_p.first);
			auto e_var = variables[e];
			nw_constr+=(-1)*e_var;
		}
		m_model.add(IloRange{m_env, 0, nw_constr, 0});

		//east-west
		IloExpr ew_constr{m_env};
		for(auto e_p=boost::out_edges(east_map[v], graph); e_p.first!=e_p.second; ++e_p.first){
			auto e = *(e_p.first);
			auto e_var = variables[e];
			ew_constr+=1*e_var;
			visited_constr+=1*e_var;
		}
		for(auto e_p=boost::out_edges(west_map[v], graph); e_p.first!=e_p.second; ++e_p.first){
			auto e = *(e_p.first);
			auto e_var = variables[e];
			ew_constr+=(-1)*e_var;
		}
		m_model.add(IloRange{m_env, 0, ew_constr, 0});

		m_model.add(IloRange{m_env, 1, visited_constr, 1});
	}

	m_cplex = IloCplex(m_model);
	//if(m_DIST_COST==0) m_cplex.setParam(IloCplex::ObjDif, 0.95*m_TURN_COST);
	std::cout << "Done." <<std::endl;
}


size_t subset::form2::Solver::eliminate_subtours(subset::SeparationMethod sm)
{
	if(sm==subset::SeparationMethod::ADVANCED) throw std::invalid_argument{"Formulation 2 does not have SeparationMethod::ADVANCED"};
	typedef boost::property<boost::vertex_name_t, grid_graph::vertex_descriptor> map_back_property;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, map_back_property, boost::no_property> comp_graph;
	std::map<grid_graph::vertex_descriptor, comp_graph::vertex_descriptor > gg_to_cg;

	comp_graph cg;
	size_t cuts_added = 0;

	//add vertices
	auto map_back = boost::get(boost::vertex_name, cg);
	for(auto v_gg: subset_vertices){
		auto v_cg = boost::add_vertex(cg);
		boost::put(map_back, v_cg, v_gg);
		gg_to_cg[v_gg] = v_cg;
	}

	//add edges
	std::map<internal_graph::edge_descriptor, IloNum> solution;
	for(auto ev: variables){
		solution[ev.first] = m_cplex.getValue(ev.second);
	}
	for(auto val_it: solution){
		if(val_it.second>0.1){
			auto v_gg = internal_to_grid_graph[boost::target(val_it.first, graph)];
			auto w_gg = internal_to_grid_graph[boost::source(val_it.first, graph)];
			if(v_gg!=w_gg){
				boost::add_edge(gg_to_cg[v_gg], gg_to_cg[w_gg], cg);
			}
		}
	}

	//calculate components
	std::vector<int> component(boost::num_vertices(cg));
	int num = boost::connected_components(cg, &component[0]);
	if(num>1){
		//create map with all components
		std::map<int, std::vector<grid_graph::vertex_descriptor>> comps;
		for(auto v_it=boost::vertices(cg); v_it.first!=v_it.second; ++v_it.first){
			auto v = *v_it.first;
			auto cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v)];
			assert(cmp_idx<=num);
			if(!comps.count(cmp_idx)) comps[cmp_idx] = std::vector<grid_graph::vertex_descriptor>();
			comps[cmp_idx].push_back(boost::get(boost::get(boost::vertex_name, cg), v));
		}

		//add constraint for each component
		for(auto comp: comps){
			std::cout << "Comp "<<comp.first<<" with "<<comp.second.size() << " elements"<<std::endl;
			IloExpr constr{m_env};
			for(auto v: comp.second){
				for(auto* v_int_map: {&north_map, &east_map, &south_map, &west_map}){
					auto v_int = (*v_int_map)[v];
					for(auto e_it=boost::out_edges(v_int, graph); e_it.first!=e_it.second; ++e_it.first){
						auto v1 = internal_to_grid_graph[boost::target(*e_it.first, graph)];
						auto v2 = internal_to_grid_graph[boost::source(*e_it.first, graph)];
						auto v1_aux = gg_to_cg[v1];
						auto v2_aux = gg_to_cg[v2];
						auto v1_cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v1_aux)];
						auto v2_cmp_idx = component[boost::get(boost::get(boost::vertex_index, cg), v2_aux)];
						if(v1_cmp_idx!=v2_cmp_idx){
							assert(solution[*e_it.first]<0.1);
							constr+=1*variables[*e_it.first];
						}
					}
				}
			}
			m_model.add(IloRange(m_env, 2, constr, IloInfinity));
			++cuts_added;
		}
		std::cout << "Solution consists of "<<num<<" components. Added "<<cuts_added<<" cuts" <<std::endl;
	} else {
		std::cout << "Solution is connected. No cuts added." <<std::endl;
	}
	return cuts_added;
}
