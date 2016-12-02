//
// Third IP formulation for subset coverage cycle cover and tour in grid graphs with distance and turn costs.
// Please note that this code was only a small part of my master's thesis and the priority was on getting working code quickly.
// The main work is done by CPLEX, hence the code does also not need to be efficient.
// 
#ifndef TURNCOSTTOURSOPTIMIZATION_FORMULATION_3_H
#define TURNCOSTTOURSOPTIMIZATION_FORMULATION_3_H

#include "../cplex.hpp"
#include "../problem_instance/grid_graph.h"
#include <ctime>
#include "./base_solver.h"

namespace subset{
	namespace form3 {
		using namespace gg;

		class Solver: public BaseSolver {


			// internal graph for IP formulation
			typedef boost::property<boost::edge_weight_t, double> FloatEdgeWeightProperty;  //integral edge weights
			typedef boost::directed_graph<boost::no_property, FloatEdgeWeightProperty> InternalGraph;
			InternalGraph internal_graph;

			//every field has a four vertex. One for each direction
			std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> north_map;
			std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> east_map;
			std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> south_map;
			std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> west_map;

			std::map<InternalGraph::vertex_descriptor, grid_graph::vertex_descriptor> internal_to_external_graph;

			std::map<InternalGraph::edge_descriptor, IloNumVar> variables;

			public:
			/**
			 * Creates the IP formulation for the instance
			 */
			Solver(grid_graph &g, std::set<coord>& subset, double turn_cost=1.0, double dist_cost=0.0, std::time_t deadline=0);


			/**
			 * Checks if the current solution is a tour and adds subtour elimination constraints if not.
			 * The return value equals the number of cuts added.
			 */
			size_t eliminate_subtours(SeparationMethod sm=SeparationMethod::ADVANCED);

			size_t cut_1(std::map<InternalGraph::edge_descriptor, IloNum>& edge_variable_assignments);
			size_t cut_2(std::map<InternalGraph::edge_descriptor, IloNum>& edge_variable_assignments);

		};

	}
}
#endif //TURNCOSTTOURSOPTIMIZATION_FORMULATION_3_H
