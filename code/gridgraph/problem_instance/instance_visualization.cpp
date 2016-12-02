
#include <iomanip>
#include "instance_visualization.h"
#include "simple_svg_1.0.0.hpp"

void gg::visualize_instance(const std::string &filename, gg::grid_graph &instance) {
    svg::Document simpleSvg(filename, svg::Layout(svg::Dimensions(10*(instance.m_max_x-instance.m_min_x+1), 10*(instance.m_max_y-instance.m_min_y+1))));
    for(auto v_p=boost::vertices(instance.grid_graph_boost); v_p.first!=v_p.second; ++v_p.first){
        auto v = *(v_p.first);
        auto v_coord = instance.getCoord(v);

        //draw
        auto stroke = svg::Stroke(0.5, svg::Color::Black);
        auto color_gray = svg::Color(169,169,169);
        auto fill = svg::Fill(color_gray);
        auto shift_x = -instance.m_min_x;
        auto shift_y = -instance.m_min_y;
        int sq = 10; //square size
        simpleSvg << svg::Rectangle(svg::Point((v_coord.first+shift_x)*sq, (v_coord.second+shift_y+1)*sq), sq, sq, fill, stroke);
    }
    simpleSvg.save();
}

void gg::visualize_instance(const std::string &filename, gg::grid_graph &instance, std::set<coord> &subset) {
    svg::Document simpleSvg(filename, svg::Layout(svg::Dimensions(10*(instance.m_max_x-instance.m_min_x+1), 10*(instance.m_max_y-instance.m_min_y+1))));
    for(auto v_p=boost::vertices(instance.grid_graph_boost); v_p.first!=v_p.second; ++v_p.first){
        auto v = *(v_p.first);
        auto v_coord = instance.getCoord(v);

        //draw
        auto stroke = svg::Stroke(0.5, svg::Color::Black);
        auto color_gray = svg::Color(169,169,169);
        auto fill = svg::Fill((subset.count(v_coord)?svg::Color::Red:color_gray));
        auto shift_x = -instance.m_min_x;
        auto shift_y = -instance.m_min_y;
        int sq = 10; //square size
        simpleSvg << svg::Rectangle(svg::Point((v_coord.first+shift_x)*sq, (v_coord.second+shift_y+1)*sq), sq, sq, fill, stroke);
    }
    simpleSvg.save();
}

void gg::visualize_instance(const std::string &filename, gg::grid_graph &instance, std::map<coord, double>& penalty) {
    svg::Document simpleSvg(filename, svg::Layout(svg::Dimensions(10*(instance.m_max_x-instance.m_min_x+1), 10*(instance.m_max_y-instance.m_min_y+1))));
    for(auto v_p=boost::vertices(instance.grid_graph_boost); v_p.first!=v_p.second; ++v_p.first){
        auto v = *(v_p.first);
        auto v_coord = instance.getCoord(v);

        //draw
        auto stroke = svg::Stroke(0.5, svg::Color::Black);
        auto color_gray = svg::Color(169,169,169);
        auto fill = svg::Fill(color_gray);
        auto shift_x = -instance.m_min_x;
        auto shift_y = -instance.m_min_y;
        int sq = 10; //square size
        simpleSvg << svg::Rectangle(svg::Point((v_coord.first+shift_x)*sq, (v_coord.second+shift_y+1)*sq), sq, sq, fill, stroke);
        //penalty label
        std::stringstream stream; stream << std::fixed << std::setprecision(2) << penalty[v_coord];
        std::string label = stream.str();
        auto text_color = svg::Color::Black;
        auto text_font = svg::Font(2.0, "Verdana");
        auto text_pos = svg::Point((v_coord.first+shift_x)*sq+sq/4, (v_coord.second+shift_y)*sq+sq/2);
        simpleSvg << svg::Text(text_pos, label,text_color, text_font);
    }
    simpleSvg.save();
}
