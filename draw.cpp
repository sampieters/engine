//
// Created by Sam Pieters on 19/02/2020.
//

#include "draw.h"

int rounddouble(double currentcolor){
    return (int)round(currentcolor);
}

img::Color convertcolor(Color currentcolor) {
    img::Color givebackcolor;
    givebackcolor.red = rounddouble(currentcolor.red * 255);
    givebackcolor.green = rounddouble(currentcolor.green * 255);
    givebackcolor.blue = rounddouble(currentcolor.blue * 255);
    return givebackcolor;
}

Line2D::Line2D(const Point2D &p1, const Point2D &p2, const Color &color) : p1(p1), p2(p2), color(color) {}

Line2D::Line2D(const Point2D &p1, const Point2D &p2, const Color &color, const double z1, const double z2)  : p1(p1), p2(p2), color(color), z1(z1), z2(z2)  {

}

std::string replace(const std::string& replacement, int iterations, const LParser::LSystem2D& l_system, std::string total){
    for (char i : replacement) {
        if(i == '+' or i == '-' or i == '(' or i == ')') {
            total.push_back(i);
        }
        else{
            total += (l_system.get_replacement(i));
        }
    }
    while (iterations != 1) {
        iterations -= 1;
        return replace(total, iterations, l_system, "");
    }
    return total;
}

using Lines2D = std::list<Line2D>;

Lines2D drawLSystem(std::string file, int size, std::vector<int> background, std::vector<double> foreground) {
    // Maak een backgroundcolor die in de .ini werd gegeven
    img::Color back(background[0], background[1], background[2]);
    Color fore(foreground[0], foreground[1], foreground[2]);

    // Zet alles klaar voor de l_parser
    LParser::LSystem2D l_system;
    std::ifstream input_stream(file);
    input_stream >> l_system;
    input_stream.close();

    // Haal alle informatie uit de.L2D file
    auto alfabet = l_system.get_alphabet();
    unsigned int iterations = l_system.get_nr_iterations();
    std::string iterator = l_system.get_initiator();
    double startangle = l_system.get_starting_angle();
    double addangle = l_system.get_angle();
    std::string bigstring;
    Point2D currentpoint{};
    Lines2D lines;

    // Grote string maken voor de lijnen te tekenen
    for (char i : iterator) {
        if(i == '+' or i == '-' or i == '(' or i== ')') {
            bigstring += i;
        } else {
            bigstring += replace(l_system.get_replacement(i), iterations-1, l_system, "");
        }
    }

    std::vector<std::pair<Point2D, double>> stack;
    // Loopen over de grote string en lijnen tekenen
    img::EasyImage image(size, size, back);
    for (char j : bigstring){
        if(j == '+') {
            startangle += addangle;
        }
        else if(j == '-') {
            startangle -= addangle;
        }
        else if(j == '(') {
            std::pair<Point2D, double> new_item;
            new_item.first = currentpoint;
            new_item.second = startangle;
            stack.push_back(new_item);
        }
        else if(j == ')') {
            currentpoint = stack.back().first;
            startangle = stack.back().second;
            stack.pop_back();
        }
        else{
            Point2D newpoint{};
            newpoint.x = currentpoint.x + std::cos(startangle * (M_PI/180));
            newpoint.y = currentpoint.y + std::sin(startangle * (M_PI/180));
            if(l_system.draw(j)){
                Line2D line = Line2D(currentpoint, newpoint, fore);
                lines.push_back(line);
            }
            currentpoint = newpoint;
        }
    }
    return lines;
}

img::EasyImage draw2DLines(Lines2D &lines, const int size, const img::Color& backgroundcolor, std::string type) {
    // (1) Cast de afbeelding naar de juiste plaats
    // 1) Bepaal Xmin, Xmax, Ymin, Ymax
    double xmin = lines.begin()->p1.x;
    double xmax = lines.begin()->p1.x;
    double ymin = lines.begin()->p1.y;
    double ymax = lines.begin()->p1.y;
    for (Line2D line:lines) {
        if(line.p1.x < xmin) {xmin = line.p1.x;}
        if(line.p2.x < xmin) {xmin = line.p2.x;}
        if(line.p1.x > xmax) {xmax = line.p1.x;}
        if(line.p2.x > xmax) {xmax = line.p2.x;}
        if(line.p1.y < ymin) {ymin = line.p1.y;}
        if(line.p2.y < ymin) {ymin = line.p2.y;}
        if(line.p1.y > ymax) {ymax = line.p1.y;}
        if(line.p2.y > ymax) {ymax = line.p2.y;}
    }
    // 2) Bereken de grootte van de Image
    // 2.1) Bepaal de Xrange
    double max = 0;
    if (xmax-xmin < ymax-ymin) {
        max = ymax-ymin;
    } else {
        max = xmax-xmin;
    }
    // 2.2) Bereken imagex en imagey
    double imagex = size * ((xmax-xmin)/max);
    double imagey = size * ((ymax-ymin)/max);

    // 3) Schaal de lijntekening
    // 3.1) Bereken de schaalfactor d
    double d = 0.95 * (imagex/(xmax-xmin));

    // 3.2) Vermenigvuldig de coördinaten van alle punten met d

    for (Line2D &line:lines) {
        line.p1.x *= d;
        line.p1.y *= d;
        line.p2.x *= d;
        line.p2.y *= d;
    }
    // 4) Verschuif de lijntekening
    // 4.1) Bereken DCx, DCy, dx, dy
    double DCx = d* ((xmin + xmax)/2);
    double DCy = d* ((ymin + ymax)/2);
    double dx = (imagex/2)-DCx;
    double dy = (imagey/2)-DCy;

    // 4.2) Tel bij alle punten van de lijntekening (dx, dy) op
    for (Line2D &line:lines) {
        line.p1.x += dx;
        line.p1.y += dy;
        line.p2.x += dx;
        line.p2.y += dy;
    }
    // (2) Beeld een willekeurige lijst met lijnen af op een image
    img::EasyImage image(imagex, imagey, backgroundcolor);
    ZBuffer ZBuffer(imagex, imagey);
    for (Line2D line:lines) {
        auto a = convertcolor(line.color);
        if(type == "ZBufferedWireframe") {
            ZBuffer.draw_zbuf_line(image, line.p1.x, line.p1.y, line.z1, line.p2.x, line.p2.y, line.z2, line.color);
        } else {
            image.draw_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y, a);
        }
    }
    return image;
}
