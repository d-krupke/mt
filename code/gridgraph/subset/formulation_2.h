//
// Second IP formulation for subset coverage cycle cover and tour in grid graphs with distance and turn costs.
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
#include <ctime>
#include "./base_solver.h"

namespace subset{
	namespace form2 {
		using namespace gg;

		class Solver: public BaseSolver {
			public:

				//internal graph representation
				typedef boost::property<boost::edge_weight_t, double> internal_graph_edge_weight_property;
				typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, internal_graph_edge_weight_property> internal_graph;
				internal_graph graph;


				//every field has 4 vertices in the internal graph, one for each direction
				std::map<grid_graph::vertex_descriptor, internal_graph::vertex_descriptor> north_map;
				std::map<grid_graph::vertex_descriptor, internal_graph::vertex_descriptor> east_map;
				std::map<grid_graph::vertex_descriptor, internal_graph::vertex_descriptor> south_map;
				std::map<grid_graph::vertex_descriptor, internal_graph::vertex_descriptor> west_map;

				//Especially used to check if two vertices of the internal graph belong to the same field
				std::map<internal_graph::vertex_descriptor, grid_graph::vertex_descriptor> internal_to_grid_graph;

				std::map<internal_graph::edge_descriptor, IloBoolVar> variables;
				std::vector<grid_graph::vertex_descriptor> subset_vertices;

			public:

				/**
				 * Creates the IP formulation for the instance
				 */
				Solver(grid_graph &instance, std::set<coord>& subset, double turn_cost=1, double dist_cost=0, std::time_t deadline=0, bool cycle_cover=true);


				/**
				 * Checks if the current solution is a tour and adds subtour elimination constraints if not.
				 * The return value equals the number of cuts added.
				 */
				size_t eliminate_subtours(SeparationMethod sm=SeparationMethod::BASIC) override;


		};

	}
}
#endif //TURNCOSTTOURSOPTIMIZATION_FORMULATION_2_H
