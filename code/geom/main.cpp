#include <iostream>
#include <set>
#include <random>
#include "auxiliary.h"
#include "ip/OptimalSolver.h"
#include "aa/ApproximateSolver.h"
#include "visualization/Visualization.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <chrono>


std::default_random_engine generator;
std::set<coord> random_instance(size_t n)
{
    std::set<coord> ret_set;
    std::uniform_real_distribution<float> distribution(0, 1);
    for (auto i = 0; i < n; ++i) {
        auto x = distribution(generator);
        auto y = distribution(generator);
        ret_set.insert({x, y});
    }
    return ret_set;
}

void write_instance_to_file(std::string filename, std::set<coord> instance)
{
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::trunc);
    if (file.is_open()) {
        for (auto p: instance) {
            file << p.first << " " << p.second << std::endl;
        }
        file.close();
    } else {
        std::cerr << "Could not write grid graph to file \"" << filename << "\"" << std::endl;
    }
}

std::set<coord> load_instance(std::string filename)
{
    std::set<coord> coords;
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line.at(0) == '#') continue;
            std::vector<std::string> tokens;
            boost::split(tokens, line, boost::is_any_of(" "));
            if (tokens.size() >= 2) {
                coords.insert(coord(std::stod(tokens[0]), std::stod(tokens[1])));
            }
        }
        file.close();
    } else {
        std::cerr << "Could not load instance from file \"" << filename << "\"" << std::endl;
    }
    return coords;
}

int main( int argc, char *argv[] )
{
    /**
    std::set<coord> instance = random_instance(100);
    ApproximateSolver apx{instance, 2, 0.1, 1};
    bool tour = true;
    OptimalSolver opt{instance, 2, !tour, 0.1, 1};
    auto apx_solution = apx.solve(tour);
    opt.solve();
    if(tour) while(opt.separate_subcycles()>0) opt.solve();
    std::cout << apx_solution.first <<std::endl;
    std::cout << opt.get_objective_value() <<std::endl;
    print_solution("/home/doms/opt_solution.svg", opt.get_solution());
    print_solution("/home/doms/aa_solution.svg", apx_solution.second);

    return 0;
     **/

    /**
    for(size_t size=550; size<=1500; size+=50){
        for(int i=0; i<10; ++i){
            std::string filename{"/home/doms/Studium/thesis-alg-2015-krupke-ma-robots/MasterThesis/Code/geom/instances/geom_"};
            filename += std::to_string(size);
            filename += "_";
            filename += std::to_string(i);
            write_instance_to_file(filename, random_instance(size));
        }
    }
    return 0;**/

    std::set<coord> instance = load_instance(argv[1]);
    size_t resolution = static_cast<size_t>(std::stoi(argv[2]));
    bool cc = argc<=5;
    auto time_begin = std::chrono::high_resolution_clock::now();
    double turn_weight = std::stod(argv[3]);
    double dist_weight = std::stod(argv[4]);
    ApproximateSolver solver{instance, resolution, dist_weight, turn_weight};
    auto solution = solver.solve(!cc);
    std::ofstream log; log.open("./results.csv", std::ios::out | std::ios::app);
    if(log.is_open()){
        log<< argv[1] << "\t" << turn_weight << "\t"<<dist_weight <<"\t"<< (!cc?"T": "CC") << '\t' << solution.first << '\t' <<  std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-time_begin).count() <<std::endl;
    } else {
        std::cerr << "Could not open logfile" <<std::endl;
        return 1;
    }
    return 0;
}