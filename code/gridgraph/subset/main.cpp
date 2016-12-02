//
// Created by doms on 5/8/16.
//

#include <stdlib.h>
#include <iostream>
#include "formulation_1.h"
#include "formulation_2.h"
#include "formulation_3.h"
#include "../problem_instance/random_instances.h"
#include <cmath>
#include <chrono>
#include <boost/program_options.hpp>
#include "../problem_instance/save_and_load.h"
#include <iostream>                      
#include <fstream>
#include "./base_solver.h"
#include <stdexcept>

namespace po = boost::program_options;

using namespace gg;
using namespace subset;

void log_result(std::string filename, std::string instance, std::string configuration, double turn_cost, double dist_cost, double result, double time, double time_creating, double time_initial_solving, double time_separation_creation, double time_separation_solving, int solve_iterations, int cuts_added, bool solved_in_time=true)
{
	std::ofstream file; file.open(filename, std::ios::app);
	file << instance << '\t' << configuration << '\t' << result << '\t' << time << '\t' << time_creating << '\t' << time_initial_solving << '\t' << time_separation_creation << '\t' << time_separation_solving << '\t' << solve_iterations << '\t' << cuts_added << '\t' << (solved_in_time?'1':'0')<< std::endl;
}

int main(int argc, char** argv) 
{
	int deadline;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("formulation,f", po::value<int>()->default_value(1), "Using formulation 1,2, or 3.")
		("tour,t", po::value<std::string>()->default_value("none"), "Tour separation: none (cycle cover), basic, advanced")
		("instance,i", po::value<std::string>(), "Instance file path")
		("deadline,d", po::value<int>(&deadline)->default_value(30*60), "Deadline in seconds")
		("logfile,l", po::value<std::string>()->default_value("./experiment_log.log"), "Path to logfile")
		("tcost", po::value<double>()->default_value(1.0), "Turn Cost")
		("dcost", po::value<double>()->default_value(0.0), "Distance Cost")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm); 


	if(vm.count("help")>0){
		std::cout << desc << std::endl;
		return 0;
	} else if(vm.count("instance")>0){
		double dist_cost = vm["dcost"].as<double>();
		double turn_cost = vm["tcost"].as<double>();
		int formulation = vm["formulation"].as<int>();
		std::string tour_separation = vm["tour"].as<std::string>();
		std::set<coord> subset;
		grid_graph g = load_instance(vm["instance"].as<std::string>(), &subset);
		std::cout << "Loaded instance with subset.size()="<<subset.size()<<std::endl;
		bool solved_in_time = true;
		auto create_begin = std::chrono::high_resolution_clock::now();

		SeparationMethod sm;
		if(vm["tour"].as<std::string>()=="none"){
			sm = SeparationMethod::NONE;
		} else if(vm["tour"].as<std::string>()=="basic"){
			sm = SeparationMethod::BASIC;
		} else {
			if(vm["tour"].as<std::string>()!="advanced") throw std::invalid_argument{"Illegal argument for --tour"};
			sm = SeparationMethod::ADVANCED;
		}
		BaseSolver* solver;
		switch(formulation){
			case 1: solver = new form1::Solver(g, subset, turn_cost, dist_cost, std::time(nullptr)+deadline); break;
			case 2: solver = new form2::Solver(g, subset, turn_cost, dist_cost, std::time(nullptr)+deadline, sm== SeparationMethod::NONE); break;
			case 3: solver = new form3::Solver(g, subset, turn_cost, dist_cost, std::time(nullptr)+deadline); break;
		}
		auto create_end = std::chrono::high_resolution_clock::now();
		auto time_creation = std::chrono::duration_cast<std::chrono::milliseconds>(create_end-create_begin);
		auto solve_cc_begin = std::chrono::high_resolution_clock::now();
		try{ solver->solve(); } catch(...){ std::cerr<< "Out of time\n";  solved_in_time=false; }
		auto solve_cc_end = std::chrono::high_resolution_clock::now();
		auto solve_cc_time = std::chrono::duration_cast<std::chrono::milliseconds>(solve_cc_end-solve_cc_begin);
		auto rounds=1;
		auto constraints_added = 0;
		if(sm == subset::SeparationMethod::NONE) {
			log_result(vm["logfile"].as<std::string>(), vm["instance"].as<std::string>(), std::string("form")+std::to_string(formulation)+"-cc", turn_cost, dist_cost,   solver->get_objective_value(), (time_creation+solve_cc_time).count(), time_creation.count(), solve_cc_time.count(), 0, 0, rounds, constraints_added);
		} else {
			std::chrono::milliseconds time_calculating_constraints{};
			std::chrono::milliseconds time_solving_constraints{};
			while(solved_in_time){
				auto calc_s_begin = std::chrono::high_resolution_clock::now();
				auto cuts_added = solver->eliminate_subtours(sm);
				auto calc_s_end = std::chrono::high_resolution_clock::now();
				time_calculating_constraints+= std::chrono::duration_cast<std::chrono::milliseconds>(calc_s_end-calc_s_begin);
				if(cuts_added==0) break;
				auto solve_s_begin = std::chrono::high_resolution_clock::now();
				try{ solver->solve(); } catch(...){ std::cerr<< "Out of time\n";  solved_in_time=false; }
				auto solve_s_end = std::chrono::high_resolution_clock::now();
				time_solving_constraints+= std::chrono::duration_cast<std::chrono::milliseconds>(solve_s_end-solve_s_begin);
				++rounds;	
				constraints_added+=cuts_added;
			}
			log_result(vm["logfile"].as<std::string>(), vm["instance"].as<std::string>(), std::string("form")+std::to_string(formulation)+"_"+vm["tour"].as<std::string>(), turn_cost, dist_cost,  solver->get_objective_value(), (time_creation+solve_cc_time+time_calculating_constraints+time_solving_constraints).count(), time_creation.count(), solve_cc_time.count(), time_calculating_constraints.count(), time_solving_constraints.count(), rounds, constraints_added);
		}
		delete solver;
	} else {				
		std::cerr<<"Invalid Instance: No instance file provided" <<std::endl;
		std::cerr<<desc<<std::endl;
		return 1;
	}


	return 0;

}
