//
// Second IP formulation for penalty coverage cycle cover and tour in grid graphs with distance and turn costs.
// Please note that this code was only a small part of my master's thesis and the priority was on getting working code quickly.
// The main work is done by CPLEX, hence the code does also not need to be efficient.
//

#ifndef TURNCOSTTOURSOPTIMIZATION_FORMULATION_2_H
#define TURNCOSTTOURSOPTIMIZATION_FORMULATION_2_H

#include <boost/graph/properties.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/named_graph.hpp>
#include "../problem_instance/grid_graph.h"
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include "../cplex.hpp"
#include "boost/multi_array.hpp"
#include "../problem_instance/distance_calculator.h"
#include "base_solver.h"

#include <ctime>
namespace penalty{
	namespace form2 {
		using namespace gg;

		class Solver: public BaseSolver {
			public:

				//internal graph representation
				typedef boost::property<boost::edge_weight_t, double> InternalGraphEdgeWeightProperty;
				typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, InternalGraphEdgeWeightProperty> InternalGraph;
				InternalGraph graph;


				//every field has 4 vertices in the internal graph, one for each direction
				std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> north_map;
				std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> east_map;
				std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> south_map;
				std::map<grid_graph::vertex_descriptor, InternalGraph::vertex_descriptor> west_map;

				//Especially used to check if two vertices of the internal graph belong to the same field
				std::map<InternalGraph::vertex_descriptor, grid_graph::vertex_descriptor> internal_to_grid_graph;

				//CPLEX
				std::map<InternalGraph::edge_descriptor, IloBoolVar> variables;

			public:

				/**
				 * Creates the IP formulation for the instance
				 */
				Solver(grid_graph &instance, std::map<coord, double>& penalties, double turn_cost=1, double dist_cost=0, std::time_t deadline=0);


				/**
				 * Checks if the current solution is a tour and adds subtour elimination constraints if not.
				 * The return value equals the number of cuts added.
				 */
				size_t eliminate_subtours(SeparationMethod sm=SeparationMethod::BASIC);



		};

	}
}
#endif //TURNCOSTTOURSOPTIMIZATION_FORMULATION_2_H
