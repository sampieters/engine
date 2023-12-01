//
// Created by Sam Pieters on 19/02/2020.
//

#ifndef ENGINE_DRAW_H
#define ENGINE_DRAW_H

#include "easy_image.h"
#include "cmath"
#include "l_parser.h"
#include "vector3d.h"
#include "Light.h"
#include "Color.h"
#include "cassert"

#include <list>
#include <fstream>
#include <limits>

class Point2D {
public:
    double x;
    double y;
};

class Line2D {
public:
    Line2D(const Point2D &p1, const Point2D &p2, const Color &color);
    Line2D(const Point2D &p1, const Point2D &p2, const Color &color, const double z1, const double z2);
    Point2D p1;
    Point2D p2;
    Color color;
    double z1;
    double z2;
};

using Lines2D = std::list<Line2D>;
int rounddouble(double currentcolor);
img::Color convertcolor(Color currentcolor);
std::string replace(const std::string& replacement, int iterations, const LParser::LSystem2D& l_system, std::string total);
Lines2D drawLSystem(std::string file, int size, std::vector<int> background, std::vector<double> foreground);
img::EasyImage draw2DLines(Lines2D &lines, const int size, const img::Color& backgroundcolor, std::string type);


#endif //ENGINE_DRAW_H
