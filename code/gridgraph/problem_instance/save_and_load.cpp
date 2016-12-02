
#include "save_and_load.h"

void gg::save_instance(grid_graph& instance, std::string filename){
	std::ofstream file;
	file.open(filename, std::ios::out|std::ios::trunc);
	if(file.is_open()) {
		for (auto v_it = boost::vertices(instance.grid_graph_boost); v_it.first != v_it.second; ++v_it.first) {
			auto v = *v_it.first;
			auto v_coord = instance.getCoord(v);
			file << v_coord.first << " " << v_coord.second << std::endl;
		}
		file.close();
	} else {
		std::cerr << "Could not write grid graph to file \"" << filename << "\"" << std::endl;
	}
}

gg::grid_graph gg::load_instance(std::string filename){
	std::set<coord> coords;
	std::ifstream file(filename);
	if(file.is_open()){
		std::string line;
		while(std::getline(file, line)){
			if(line.empty()) continue;
			if(line.at(0)=='#') continue;
			std::vector<std::string> tokens;
			boost::split(tokens, line, boost::is_any_of(" "));
			if(tokens.size()>=2){
				coords.insert(std::pair<int, int>(std::stoi(tokens[0]), std::stoi(tokens[1])));
			}
		}
		file.close();
	} else {
		std::cerr << "Could not load grid graph from file \"" << filename << "\"" << std::endl;
	}
	return grid_graph(coords);
}

void gg::save_instance(gg::grid_graph& instance, std::set<coord>& subset, std::string filename){
	std::ofstream file;
	file.open(filename, std::ios::out|std::ios::trunc);
	if(file.is_open()) {
		for (auto v_it = boost::vertices(instance.grid_graph_boost); v_it.first != v_it.second; ++v_it.first) {
			auto v = *v_it.first;
			auto v_coord = instance.getCoord(v);
			file << v_coord.first << " " << v_coord.second << " " << subset.count(v_coord) << std::endl;
		}
		file.close();
	} else {
		std::cerr << "Could not write grid graph to file \"" << filename << "\"" << std::endl;
	}
}

gg::grid_graph gg::load_instance(std::string filename, std::map<coord, double> &penalty) {
	std::set<coord> coords;
	std::ifstream file(filename);
	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {
			if (line.empty()) continue;
			if (line.at(0) == '#') continue;
			std::vector<std::string> tokens;
			boost::split(tokens, line, boost::is_any_of(" "));
			if (tokens.size() >= 3) {
				coord c(std::stoi(tokens[0]), std::stoi(tokens[1]));
				coords.insert(c);
				penalty[c] = std::stod(tokens[2]);
			} else {
				assert(false);
			}
		}
		file.close();
	} else {
		std::cerr << "Could not load grid graph from file \"" << filename << "\"" << std::endl;
	}
	return grid_graph(coords);
}

void gg::save_instance(gg::grid_graph &instance, std::map<coord, double> &penalty, std::string filename) {
	std::ofstream file;
	file.open(filename, std::ios::out | std::ios::trunc);
	if (file.is_open()) {
		for (auto v_it = boost::vertices(instance.grid_graph_boost); v_it.first != v_it.second; ++v_it.first) {
			auto v = *v_it.first;
			auto v_coord = instance.getCoord(v);
			file << v_coord.first << " " << v_coord.second << " " << penalty[v_coord] << std::endl;
		}
		file.close();
	} else {
		std::cerr << "Could not write grid graph to file \"" << filename << "\"" << std::endl;
	}
}

gg::grid_graph gg::load_instance(std::string filename, std::set<coord> *subset) {
	std::set<coord> coords;
	std::ifstream file(filename);
	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {
			if (line.empty()) continue;
			if (line.at(0) == '#') continue;
			std::vector<std::string> tokens;
			boost::split(tokens, line, boost::is_any_of(" "));
			if (tokens.size() >= 3) {
				coord c(std::stoi(tokens[0]), std::stoi(tokens[1]));
				coords.insert(c);
				auto in_subset = std::stoi(tokens[2]);
				if(in_subset!=0 && in_subset!=1){
					std::cerr << "File is wrong coded. Third column is only allowed to consist of '0' and '1'" << std::endl;
				}
				if (in_subset==1) {
					subset->insert(c);
				}
			} else {
				assert(false);
			}
		}
		file.close();
	} else {
		std::cerr << "Could not load grid graph from file \"" << filename << "\"" << std::endl;
	}
	return grid_graph(coords);
}
