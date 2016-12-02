#include "./grid_graph.h"
#include "./random_instances.h"
#include "./instance_visualization.h"
#include "./save_and_load.h"
#include <cmath>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace gg;

grid_graph create_dense_instance(size_t size){
	size_t bounding_box = std::round(std::sqrt(3*size));
	return create_random_grid_graph(size, bounding_box, bounding_box, 0.001, 5, 4);
}

grid_graph create_sparse_instance(size_t size){
	size_t bounding_box = std::round(std::sqrt(3*size));
	return create_random_grid_graph(size, bounding_box, bounding_box, 0.01, 3, 2);
}



int main(int argc, char** argv){

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("problem,p", po::value<std::string>()->default_value("fc"), "Create an instance for problem 'fc', 'subset', or 'penalty'")
		("size,s", po::value<unsigned int>(), "Size of instance")
		("load,l", po::value<std::string>(), "Load instance")
		("dense", "Create a dense instance")
		("sparse", "Create a sparse instance")
		("subsetsize", po::value<unsigned int>(), "Size of subset")
		("avg", po::value<double>(), "Average")
		("dev", po::value<double>(), "Deviation")
		("output,o", po::value<std::string>(), "Output file")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm); 

	if(vm.count("help")>0){
		std::cout << desc << std::endl;
		return 0;
	} else {
		grid_graph gg;
		if(vm.count("load")>0){
			gg = load_instance(vm["load"].as<std::string>());
		} else if(vm.count("dense")>0){
			gg = create_dense_instance(vm["size"].as<unsigned int>());
		} else {
			gg = create_sparse_instance(vm["size"].as<unsigned int>());
		}

		if(vm["problem"].as<std::string>()=="s"){
			if(vm.count("subsetsize")==0) { return 1; }
			auto subset = create_random_subset(gg, vm["subsetsize"].as<unsigned int>());
			save_instance(gg, subset, vm["output"].as<std::string>());
			return 0;
		} else if(vm["problem"].as<std::string>()=="p") {
			if(vm.count("avg")+vm.count("dev")<2) return 1;
			auto penalty = create_random_penalty(gg, vm["avg"].as<double>(), vm["dev"].as<double>());
			save_instance(gg, penalty, vm["output"].as<std::string>());
			return 0;
		} else {
			save_instance(gg, vm["output"].as<std::string>());

			return 0;
		}
	}

	/**
	int size = std::stoi(argv[1]);
	int max_x = std::stoi(argv[2]);
	int max_y = std::stoi(argv[3]);
	double delete_prob = std::stod(argv[4]);
	double inc_direct_nbr = std::stod(argv[5]);
	double inc_diag_nbr = std::stod(argv[6]);
	grid_graph instance = create_random_grid_graph(size, max_x, max_y, delete_prob, inc_direct_nbr, inc_diag_nbr);
	//grid_graph instance = create_random_grid_graph_from_png("./png/03.png", 1000, 0.01f, 3.0f, 2.0f);
	//auto penalties = create_random_penalty_from_png(instance, "./png/_02.png", 2, 1);
	visualize_instance(std::string("test.svg"), instance);
	save_instance(instance, "test.gg");
	**/

}
