
#include "Visualization.h"
#include "./simple_svg_1.0.0.hpp"

void print_solution(std::string filename, std::vector<std::pair<coord, coord>> solution)
{
    auto coord_to_svg = [](coord p){return svg::Point{5*p.first, 5*p.second};};
    svg::Document svgfile{filename, svg::Layout{svg::Dimensions{5,5}}};

    auto stroke_edge = svg::Stroke{0.01f, svg::Color::Red};

    for(auto e: solution){
        svgfile << svg::Line(coord_to_svg(e.first), coord_to_svg(e.second), stroke_edge);
        svgfile << svg::Circle(coord_to_svg(e.first), 0.03, svg::Fill(svg::Color::Red), svg::Stroke{0.01f, svg::Color::Red});
        svgfile << svg::Circle(coord_to_svg(e.second), 0.03, svg::Fill(svg::Color::Red), svg::Stroke{0.01f, svg::Color::Red});
    }
    svgfile.save();
}
