
#ifndef GRIDGRAPH_PENALTY_BASE_SOLVER_H_H
#define GRIDGRAPH_PENALTY_BASE_SOLVER_H_H

#include "../cplex.hpp"
#include "../problem_instance/grid_graph.h"
#include <map>
#include <ctime>

namespace penalty
{
	using namespace gg;
	using gg_vd = grid_graph::boost_grid_graph_t::vertex_descriptor;
	enum class SeparationMethod { BASIC, ADVANCED, NONE };

	class BaseSolver {

		protected:
			//instance
			const double m_DIST_COST;
			const double m_TURN_COST;
			const std::time_t m_DEADLINE;

			//Problem instance
			grid_graph& m_instance;
			std::map<coord, double> m_penalties;

			//CPLEX
			IloEnv m_env;
			IloModel m_model;
			IloCplex m_cplex;
			IloObjective m_objective;

			double lowerbound = 0.0;
			bool cplex_params_set = false;
			bool error_occured = false;

		public:
			BaseSolver(grid_graph &instance, std::map<coord, double> penalties, double turn_cost=0.0, double dist_cost=0.0, std::time_t deadline=0):
				m_DIST_COST{dist_cost}, m_TURN_COST{turn_cost},  m_DEADLINE{deadline}, m_instance{instance}, m_penalties{penalties}, m_model{m_env}
			{
				std::cout << "Instance with turn_cost="<<turn_cost<<" and dist_cost="<<dist_cost<<std::endl;
			}


			IloAlgorithm::Status solve()
			{
				if(error_occured){
					std::cerr << "Calling solve with error alread occured" <<std::endl;
				}
				if(!cplex_params_set){
					m_cplex.setParam(IloCplex::EpInt, 0);
					m_cplex.setParam(IloCplex::EpGap, 0);
					m_cplex.setParam(IloCplex::EpOpt, 1e-9);
					m_cplex.setParam(IloCplex::EpAGap, 0);
					//m_cplex.setParam(IloCplex::NumericalEmphasis, 1);
					cplex_params_set = true;
				}
				if(m_DEADLINE!=0){
					std::time_t now = std::time(nullptr);
					long resttime = m_DEADLINE-now;
					if(resttime<=0){
						std::cout << "Abort solve because no time left." <<std::endl;
						throw (-1);
					}
					m_cplex.setParam(IloCplex::TiLim, resttime);
				}
				std::cout << "Solving IP...." <<std::endl;
				//Solve
				try{
					if (!m_cplex.solve()) {
						error_occured = true;
						m_env.error() << "Failed to optimize LP: " << m_cplex.getStatus() << std::endl;
						throw (-1);
					}
				} catch (IloException e){
					std::cerr << "Error while solving: "<<e <<std::endl;
				}
				std::cout << "Solved IP to optimality with solution: "<<m_cplex.getObjValue()<<std::endl;
				if(m_DEADLINE<std::time(nullptr)) throw (-1);
				lowerbound = m_cplex.getObjValue();
				return m_cplex.getStatus();
			}


			double get_objective_value(){
				if(error_occured) return -2;
				constexpr double RETURN_VALUE_IF_OUT_OF_TIME = -1;
				if(m_DEADLINE<std::time(nullptr)) return RETURN_VALUE_IF_OUT_OF_TIME;
				double ret;
				try{
					ret = m_cplex.getObjValue();
				} catch (...){
					std::cerr << "ERROR: Exception thrown while querying objective value. Return 0." << std::endl;
					return RETURN_VALUE_IF_OUT_OF_TIME;
				}
				return ret;
			}

			/**
			 * Prohibits the subtours in the current solution. Use it until it returns zero to obtain a tour.
			 * The returned value returns the amount of added constraints.
			 * There are at most two different separation methods: Basic and Advanced where Advanced has an additional
			 * constraint for formulation 1 and 3. Formulation 2 only has Basic.
			 **/
			virtual size_t eliminate_subtours(SeparationMethod sm)=0;

			virtual ~BaseSolver()=default;
	};
}
#endif //GRIDGRAPH_PENALTY_BASE_SOLVER_H_H
