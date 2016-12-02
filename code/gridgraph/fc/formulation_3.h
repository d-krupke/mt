//
// Third IP formulation for full coverage cycle cover and tour in grid graphs with distance and turn costs.
// Please note that this code was only a small part of my master's thesis and the priority was on getting working code quickly.
// The main work is done by CPLEX, hence the code does also not need to be efficient.
// 

#ifndef GG_FC_F3_H
#define GG_FC_F3_H

#include "../cplex.hpp"
#include "../problem_instance/grid_graph.h"
#include <ctime>
#include "./base_solver.h"

namespace fc{
	namespace form3 {
		using namespace gg;

		class Solver: public BaseSolver {

			// internal graph for IP formulation
			typedef boost::property<boost::edge_weight_t, float> EdgeWeightProperty_int;  //integral edge weights
			typedef boost::directed_graph<boost::no_property, EdgeWeightProperty_int> aux_graph;
			aux_graph internal_graph;

			//every field has a four vertex. One for each direction
			std::map<grid_graph::vertex_descriptor, aux_graph::vertex_descriptor> north_map;
			std::map<grid_graph::vertex_descriptor, aux_graph::vertex_descriptor> east_map;
			std::map<grid_graph::vertex_descriptor, aux_graph::vertex_descriptor> south_map;
			std::map<grid_graph::vertex_descriptor, aux_graph::vertex_descriptor> west_map;

			std::map<aux_graph::vertex_descriptor, grid_graph::vertex_descriptor> internal_to_external_graph;

			std::map<aux_graph::edge_descriptor, IloNumVar> variables;


			public:
			/**
			 * Creates the IP formulation for the instance
			 */
			Solver(grid_graph &g, double turn_cost=1, double dist_cost=0, std::time_t deadline=0);


			/**
			 * Checks if the current solution is a tour and adds subtour elimination constraints if not.
			 * The return value equals the number of cuts added.
			 */
			size_t eliminate_subtours(SeparationMethod sm=SeparationMethod::ADVANCED);

			size_t cut_1(std::map<aux_graph::edge_descriptor, IloNum>& solution);
			size_t cut_2(std::map<aux_graph::edge_descriptor, IloNum>& solution);

		};

	}
}
#endif //GG_FC_F3_H
