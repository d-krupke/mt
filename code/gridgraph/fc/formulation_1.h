//
// First IP formulation for full coverage cycle cover and tour in grid graphs with distance and turn costs.
// Please note that this code was only a small part of my master's thesis and the priority was on getting working code quickly.
// The main work is done by CPLEX, hence the code does also not need to be efficient.
// 
#ifndef GG_FC_F1_H
#define GG_FC_F1_H

#include "../cplex.hpp"
#include "../problem_instance/grid_graph.h"
#include <ctime>
#include "./base_solver.h"

namespace fc{
	namespace form1 {
		using namespace gg;

		//a variable is identified by the two fields it connects and the intermediate field
		struct var_index: std::tuple<gg_vd, gg_vd, gg_vd>{
			var_index():std::tuple<gg_vd, gg_vd, gg_vd>{}{}
			var_index(gg_vd u, gg_vd v, gg_vd w): std::tuple<gg_vd, gg_vd, gg_vd>{std::min(u,w),v,std::max(u,w)}{}
		};

		class Solver: public BaseSolver {
			std::map<var_index, IloIntVar> variable_map; //Saves the corresponding variable for each field passing

			public:
			// Creates the base IP
			Solver(grid_graph& instance, double turn_cost=1.0, double dist_cost=0.0, std::time_t deadline=0);

            //Adds subcycle elimination constraints to obtain a tour. Has to be called until it returns 0. ADVANCED uses both cuts, BASIC only the second one.
			size_t eliminate_subtours(SeparationMethod sm=SeparationMethod::ADVANCED);

			// A simple but insufficient subcycle elimination constraint to obtain tours
			size_t cut_1(std::map<var_index, IloNum>& edge_variable_assignments);
            // An advanced and sufficient subcycle elimination constraint to obtain tours. Not always as efficient as the simple one.
			size_t cut_2(std::map<var_index, IloNum>& edge_variable_assignments);

		};

	}
}
#endif // GG_FC_F1_H
