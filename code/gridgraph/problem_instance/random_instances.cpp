#include "random_instances.h"
#include <chrono>

gg::grid_graph gg::create_random_grid_graph_from_png(std::string filename, size_t desired_size, float thinout_fac, float inc_direct_nbr, float inc_diag_nbr) {
    std::vector<coord> coords;
    std::set<coord> coord_set;
/**
    png::image<png::rgb_pixel> img(filename);


    const auto height = img.get_height();
    const auto width = img.get_width();

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            coord_set.insert(coord{x, y});
        }
    }

    unsigned char red;
    unsigned char green;
    unsigned char blue;

    // 1: Full neighborhood, <<1: sparse neighborhood
    auto neighbor_value = [&](coord c) -> float {
        float val = 1;
        auto x = c.first;
        auto y = c.second;
        auto nbrs = {coord(x, y - 1), coord(x - 1, y), coord(x + 1, y), coord(x, y + 1)};
        auto diag_nbrs = {coord{x + 1, y + 1}, coord{x - 1, y + 1}, coord{x + 1, y - 1}, coord{x - 1, y - 1}};
        for (auto n: nbrs) {
            val *= (coord_set.count(n) ? 1 : 1 + inc_direct_nbr);
        }
        for (auto n: diag_nbrs) {
            val *= (coord_set.count(n) ? 1 : 1 + inc_diag_nbr);
        }
        return val;
    };


    while (coord_set.size() > desired_size) {
        std::srand(unsigned(std::time(0)));
		coords.clear();
		for(auto c: coord_set){
			coords.push_back(c);
		}
        std::random_shuffle(coords.begin(), coords.end());

        for (auto c: coords) {
            if (!coord_set.count(c)) continue;
            auto x = c.first;
            auto y = c.second;
            auto pixel = img.get_pixel(x, y);
            red = pixel.red;
            green = pixel.green;
            blue = pixel.blue;
            float Y = 0.299f * red + 0.587f * green + 0.114f * blue;
            if ((std::rand() % 2550) * thinout_fac * neighbor_value(c) > Y * 20) {
                coord_set.erase(c);
                grid_graph gg{coord_set};
                if (!gg.is_connected()) {
                    coord_set.insert(c);
                } else {
					if(coord_set.size()%100==0){
						std::cout<<  coord_set.size() <<std::endl;
					}
					if(coord_set.size()<=desired_size) break;
				}
            }
        }
    }**/
    return grid_graph{coord_set};
}

std::set<coord> gg::create_random_subset(grid_graph &instance, size_t desired_size) {
    std::srand(unsigned(std::time(0)));
    std::set<coord> subset;

    std::vector<coord> fields;
    for (auto v_p = boost::vertices(instance.grid_graph_boost); v_p.first != v_p.second; ++v_p.first) {
        auto v = *(v_p.first);
        auto v_coord = instance.getCoord(v);
        fields.push_back(v_coord);
    }
    std::shuffle(fields.begin(), fields.end(), std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
    for (auto f: fields) {
        if (subset.size() >= desired_size) break;
        subset.insert(f);
    }
    return subset;
}

std::map<coord, double> gg::create_random_penalty(grid_graph &instance, double mean, double stddev) {
    std::default_random_engine generator{std::chrono::system_clock::now().time_since_epoch().count()};
    std::normal_distribution<double> distribution(mean, stddev);

    std::map<coord, double> ret;
    for (auto v_p = boost::vertices(instance.grid_graph_boost); v_p.first != v_p.second; ++v_p.first) {
        double penalty = std::max(0.0, std::round(100*distribution(generator))/100.0);//>0 and round to *.XX
        ret[instance.getCoord(*(v_p.first))] = penalty;
    }
    return ret;
}

gg::grid_graph gg::create_random_grid_graph(size_t desired_size, int max_x, int max_y, float delete_prob_perc, float inc_direct_nbr, float inc_diag_nbr) {
    std::cout << "Creating a random grid graph...\n";

    std::vector<coord> coords;
    std::set<coord> coord_set;
    for (int x = 0; x < max_x; ++x) {
        for (int y = 0; y < max_y; ++y) {
            coords.push_back(coord{x, y});
            coord_set.insert(coord{x, y});
        }
    }
    // 1: Full neighborhood, <<1: sparse neighborhood
    auto neighbor_value = [&](coord c) {
        float val = 1;
        auto x = c.first;
        auto y = c.second;
        auto nbrs = {coord(x, y - 1), coord(x - 1, y), coord(x + 1, y), coord(x, y + 1)};
        auto diag_nbrs = {coord{x + 1, y + 1}, coord{x - 1, y + 1}, coord{x + 1, y - 1}, coord{x - 1, y - 1}};
        for (auto n: nbrs) {
            val *= (coord_set.count(n) ? 1 : 1 + inc_direct_nbr);
        }
        for (auto n: diag_nbrs) {
            val *= (coord_set.count(n) ? 1 : 1 + inc_diag_nbr);
        }
        return val;
    };

    while (coord_set.size() > desired_size) {
        std::srand(unsigned(std::time(0)));
        std::shuffle(coords.begin(), coords.end(), std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
        for (auto c: coords) {
            if (std::rand() % 1000 < 1000 * delete_prob_perc * neighbor_value(c)) {
                // std::cout << "Delete"<<std::endl;
                coord_set.erase(c);
                grid_graph gg{coord_set};
                if (!gg.is_connected()) {
                    // std::cout << "Unconnected. Reinsert." <<std::endl;
                    coord_set.insert(c);
                } else {
					if(coord_set.size()<= desired_size) break;
				}
            }
        }
    }
    std::cout << "Found connected grid graph of size "<<coord_set.size()<<"\n";
    return grid_graph{coord_set};
}

std::set<coord> gg::create_random_subset_from_png(gg::grid_graph &instance, std::string filename, int desired_subset_size, float margin) {
    std::set<coord> subset;
	/**
    png::image<png::rgb_pixel> img(filename);
    std::vector<coord> fields;
    for (auto v_p = boost::vertices(instance.grid_graph_boost); v_p.first != v_p.second; ++v_p.first) {
        auto v = *(v_p.first);
        auto v_coord = instance.getCoord(v);
        fields.push_back(v_coord);
    }
    std::srand(unsigned(std::time(0)));
    while (subset.size() < desired_subset_size && !fields.empty()) {
        std::random_shuffle(fields.begin(), fields.end());
        std::vector<coord> next_round_fields;
        for (auto c: fields) {
            if (c.first > img.get_width() || c.second > img.get_height() || c.first < 0 || c.second < 0) continue;//TODO maybe warning

            auto pixel = img.get_pixel(static_cast<size_t>(c.first), static_cast<size_t>(c.second));
            float Y = 0.299f * pixel.red + 0.587f * pixel.green + 0.114f * pixel.blue;
            if ((std::rand() % 2550) < Y * 10 * margin) {
                subset.insert(c);
            } else {
                next_round_fields.push_back(c);
            }

        }
        fields = next_round_fields;
    }**/
    return subset;
}

std::map<coord, double> gg::create_random_penalty_from_png(gg::grid_graph &instance, std::string filename, double mean, double stddev) {
    std::map<coord, double> penalties;/**
    png::image<png::rgb_pixel> img(filename);
    std::default_random_engine generator;

    for (auto v_p = boost::vertices(instance.grid_graph_boost); v_p.first != v_p.second; ++v_p.first) {
        auto v = *(v_p.first);
        auto v_coord = instance.getCoord(v);
        if (v_coord.first > img.get_width() || v_coord.second > img.get_height() || v_coord.first < 0 || v_coord.second < 0) {
            penalties[v_coord] = 0.0;
        } else {
            auto pixel = img.get_pixel(static_cast<size_t>(v_coord.first), static_cast<size_t>(v_coord.second));
            float brightness = 0.299f * pixel.red + 0.587f * pixel.green + 0.114f * pixel.blue;//0...255
            std::normal_distribution<double> distribution(mean * (brightness/255), stddev * (brightness/255));
            penalties[v_coord] = std::max(0.0, distribution(generator));
        }
    }**/
    return penalties;
}

gg::grid_graph gg::create_random_grid_graph(size_t desired_size) {
    return create_random_grid_graph(desired_size, std::round(std::sqrt(desired_size) * 2), std::round(std::sqrt(desired_size) * 2));
}
