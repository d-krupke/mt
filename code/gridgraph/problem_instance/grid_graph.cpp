
#include <iostream>
#include "grid_graph.h"

/**
bool operator<(const coord& c1, const coord& c2) {
    //coord& c1 = *this;
    return c1.x<c2.x || (c1.x==c2.x && c1.y<c2.y);
}**/

//#include "png++/png.hpp"

namespace gg {

    grid_graph::grid_graph(std::set<coord> &coords) {
        build_graph(coords);
    }

    void grid_graph::build_graph(std::set<coord> &coords){
        if(coords.empty()){
            std::cout << "WARNING: Building empty grid graph. Probably you have inserted a bad file or something similiar." << std::endl;
        }
        auto vertex_name_map = boost::get(boost::vertex_name, grid_graph_boost);
        for (auto c = coords.begin(); c != coords.end(); ++c) {

            //calculate bounding box of grid graph
            if(c==coords.begin()){
                m_min_x=m_max_x=c->first;
                m_min_y=m_max_y=c->second;
            } else {
                m_min_x=std::min(m_min_x, c->first);
                m_max_x=std::max(m_max_x, c->first);
                m_min_y=std::min(m_min_y, c->second);
                m_max_y=std::max(m_max_y, c->second);
            }

            auto v = boost::add_vertex(grid_graph_boost);
            boost::put(vertex_name_map, v, *c );
            coord_to_vertex[*c] = v;
            //.insert(*c, v);
            //std::cout << "add vertex "<<v<<" at position "<< c->first<< " "<<c->second<< " " << getCoord(v).first << std::endl;
            coord c_north{c->first, c->second + 1};
            coord c_east{c->first + 1, c->second};
            coord c_south{c->first, c->second - 1};
            coord c_west{c->first - 1, c->second};
            if (coord_to_vertex.count(c_north) ) {
                boost::add_edge(v, getVertex(c_north), grid_graph_boost );
            }
            if (coord_to_vertex.count(c_east) ) {
                boost::add_edge(v, getVertex(c_east), grid_graph_boost );
            }
            if (coord_to_vertex.count(c_south) ) {
                boost::add_edge(v, getVertex(c_south), grid_graph_boost );
            }
            if (coord_to_vertex.count(c_west) ) {
                boost::add_edge(v, getVertex(c_west), grid_graph_boost );
            }
        }
    }

    void grid_graph::load_from_file(const std::string filename) {
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
            build_graph(coords);
        } else {
            std::cerr << "Could not load grid graph from file \"" << filename << "\"" << std::endl;
        }
    }

    void grid_graph::write_to_file(const std::string filename) {
        std::ofstream file;
        file.open(filename, std::ios::out|std::ios::trunc);
        if(file.is_open()) {
            auto coord_map = boost::get(boost::vertex_name, grid_graph_boost);
            for (auto v_it = boost::vertices(grid_graph_boost); v_it.first != v_it.second; ++v_it.first) {
                auto v = *v_it.first;
                auto v_coord = boost::get(coord_map, v);
                file << v_coord.first << " " << v_coord.second << std::endl;
            }
            file.close();
        } else {
            std::cerr << "Could not write grid graph to file \"" << filename << "\"" << std::endl;
        }
    }

    void grid_graph::load_from_bitmap(const std::string filename, int threshold) {
		/**
        std::set<coord> coords;
        png::image< png::rgb_pixel > img(filename);
        const auto height = img.get_height();
        const auto width = img.get_width();

        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                auto pixel = img.get_pixel(x,y);
                float Y = 0.299f * pixel.red + 0.587f * pixel.green + 0.114f * pixel.blue;
                if(Y>threshold){
                    coords.insert(coord{x,y});
                }
            }
        }
        build_graph(coords);**/
    }




}
