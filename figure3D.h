//
// Created by Sam Pieters on 04/03/2020.
//

#ifndef ENGINE_FIGURE3D_H
#define ENGINE_FIGURE3D_H

#include "vector3d.h"
#include "draw.h"
#include "ini_configuration.h"

#include <algorithm>
#include <vector>
#include <list>
#include <utility>
#include <fstream>
#include <iostream>

using Lines2D = std::list<Line2D>;

class Face {
public:
    //De indexen refereren naar
    //punten in de ‘points’ vector
    //van de Figure-klasse
    std::vector<int> point_indexes;
    Face();
    Face(std::vector<int> point_indexes);
};

class figure3D {
public:
    figure3D();
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    //New colors for light
    Color ambientReflection;
    Color diffuseReflection;
    Color specularReflection;
    double reflectionCoefficient;

    static figure3D make_cube();
    static figure3D make_tetrahedron();
    static figure3D make_octahedron();
    static figure3D make_icosahedron();
    static figure3D make_dodecahedron();
    static figure3D make_buckyball();
    static figure3D make_sphere(const int n);
    static figure3D make_cone(const int n, const double height);
    static figure3D make_cylinder(const int n, const double height);
    static figure3D make_torus(const double r, const double R, const int n, const int m);
    static figure3D draw3LSystem(std::string inputfile);
    static std::list<figure3D> make_menger(figure3D& fig, std::list<figure3D> fractals, const int nr_iterations, int current_it);
};

typedef std::list<figure3D> Figures3D;
Figures3D generateFractal(figure3D& fig, Figures3D& fractal, const int nr_iterations, int current_it, double scale);
img::EasyImage line_drawings(const ini::Configuration &configuration);
img::EasyImage draw3DLines(Figures3D &figure, const int size, const img::Color& backgroundcolor, Lights3D lights);
Matrix scaleFigure(const double scale);
Matrix rotateX(double angle);
Matrix rotateY(double angle);
Matrix rotateZ(double angle);
Matrix translate(const Vector3D &vector);
Matrix eyePointTrans(const Vector3D &eyepoint);
void applyTransformation(figure3D& figure, const Matrix& matrix);
void applyTransformation(Figures3D &figure, const Matrix &);
void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
Lines2D doProjection(const Figures3D &);
Point2D doProjection(const Vector3D &point, const double d);
std::vector<Face> triangulate(const Face& face);
std::string replaceD3L(const std::string& replacement, int iterations, const LParser::LSystem3D& l_system, std::string total);
void generateThickFigure(const figure3D &lineDrawing, Figures3D &resultingFigures, const double r, const int n, const int m);

#endif //ENGINE_FIGURE3D_H
