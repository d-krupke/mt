//
// First IP formulation for penalty coverage cycle cover and tour in grid graphs with distance and turn costs.
// Please note that this code was only a small part of my master's thesis and the priority was on getting working code quickly.
// The main work is done by CPLEX, hence the code does also not need to be efficient.
//
#ifndef GRIDGRAPH_PENALTY_FORMULATION_1_H
#define GRIDGRAPH_PENALTY_FORMULATION_1_H

#include "../cplex.hpp"
#include "../problem_instance/grid_graph.h"
#include "base_solver.h"
#include <map>
#include <ctime>

namespace penalty{

	namespace form1 {
		using namespace gg;

		class Solver: public BaseSolver {

			//a variable is identified by the two fields it connects and the intermediate field
			struct var_index : std::tuple<gg_vd, gg_vd, gg_vd> {
				var_index() : std::tuple<gg_vd, gg_vd, gg_vd>{} {}
				var_index(gg_vd u, gg_vd v, gg_vd w) : std::tuple<gg_vd, gg_vd, gg_vd>{std::min(u, w), v, std::max(u, w)} {}
			};

			//CPLEX
			std::map<var_index, IloIntVar> variable_map;
			std::map<gg_vd, IloBoolVar> penalty_paid_variable_map;

			//Problem instance
			std::map<gg_vd, double> m_penalties_vd;

			public:
			Solver(grid_graph& m_instance, std::map<coord, double>& penalties, double turn_cost=1.0, double dist_cost=0.0, std::time_t deadline=0);


			size_t eliminate_subtours(SeparationMethod sm=SeparationMethod::ADVANCED);

			size_t cut_1(std::map<var_index, IloNum>& edge_variable_assignments, std::map<gg_vd, IloNum>& penalty_variable_assigments);
			
			size_t cut_2(std::map<var_index, IloNum>& edge_variable_assignments, std::map<gg_vd, IloNum>& penalty_variable_assigments);

		};

	}
}


#endif //GRIDGRAPH_PENALTY_FORMULATION_1_H
