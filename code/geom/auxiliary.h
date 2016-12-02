
#ifndef GEOM_AUXILIARY_H
#define GEOM_AUXILIARY_H

#include <utility>
#include <assert.h>
#include <cmath>
#include <iostream>
# define M_PI           3.14159265358979323846  /* pi */

using coord = std::pair<float, float>;
using conf = std::pair<coord, float>;

inline double dist(coord a, coord b){
    auto x = a.first-b.first;
    auto y = a.second- b.second;
    return std::sqrt(x*x+y*y);
}

inline double normalize_to_zero_twopi(double a){
    while(a<0) a+=2*M_PI;
    while(a>2*M_PI) a-=2*M_PI;
    return a;
}

inline double angle_diff(double a, double b){
    a = normalize_to_zero_twopi(a);
    b = normalize_to_zero_twopi(b);
    auto x = std::max(a,b);
    auto y = std::min(a,b);
    auto diff =x-y;
    assert(diff>=0);
    if(diff<= M_PI) return diff;
    return std::max((2*M_PI+y)-x, 0.0);
}

inline double costs(const conf a, const conf b, const double dist_weight, const double turn_weight){
    coord ab{b.first.first-a.first.first, b.first.second-a.first.second};
    coord ba{a.first.first-b.first.first, a.first.second-b.first.second};
    auto angle_a = std::atan2(ab.second, ab.first);
    auto angle_b = std::atan2(ba.second, ba.first);

    auto val =  dist_weight*dist(a.first, b.first)+turn_weight*(angle_diff(angle_a, a.second)+angle_diff(angle_b, b.second));
    //std::cout << "("<<a.first.first<<", "<<a.first.second<<", "<<a.second<<") ("<<b.first.first<<", "<<b.first.second<<", "<<b.second<<") =" << val <<std::endl;
    return val;
}

#endif //GEOM_AUXILIARY_H
