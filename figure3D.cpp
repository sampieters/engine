//
// Created by Sam Pieters on 04/03/2020.
//

#include "figure3D.h"

figure3D figure3D::make_cube() {
    Vector3D point1 = Vector3D::point(1, -1, -1);
    Vector3D point2 = Vector3D::point(-1, 1, -1);
    Vector3D point3 = Vector3D::point(1, 1, 1);
    Vector3D point4 = Vector3D::point(-1, -1, 1);
    Vector3D point5 = Vector3D::point(1, 1, -1);
    Vector3D point6 = Vector3D::point(-1, -1, -1);
    Vector3D point7 = Vector3D::point(1, -1, 1);
    Vector3D point8 = Vector3D::point(-1, 1, 1);
    std::vector<Vector3D> cube_points = {point1, point2, point3, point4, point5, point6, point7, point8};
    Face line1({0, 4, 2, 6});
    Face line2({4, 1, 7, 2});
    Face line3({1, 5, 3, 7});
    Face line4({5, 0, 6, 3});
    Face line5({6, 2, 7, 3});
    Face line6({0, 5, 1, 4});
    std::vector<Face> cube_lines = {line1, line2, line3, line4, line5, line6};
    figure3D figure;
    figure.points = cube_points;
    figure.faces = cube_lines;
    return figure;
}

figure3D figure3D::make_tetrahedron() {
    Vector3D point1 = Vector3D::point(1, -1, -1);
    Vector3D point2 = Vector3D::point(-1, 1, -1);
    Vector3D point3 = Vector3D::point(1, 1, 1);
    Vector3D point4 = Vector3D::point(-1, -1, 1);
    std::vector<Vector3D> tetrahedron_points = {point1, point2, point3, point4};
    Face line1({0, 1, 2});
    Face line2({1, 3, 2});
    Face line3({0, 3, 1});
    Face line4({0, 2, 3});
    std::vector<Face> tetrahedron_lines = {line1, line2, line3, line4};
    figure3D figure;
    figure.points = tetrahedron_points;
    figure.faces = tetrahedron_lines;
    return figure;
}

figure3D figure3D::make_octahedron() {
    Vector3D point1 = Vector3D::point(1, 0, 0);
    Vector3D point2 = Vector3D::point(0, 1, 0);
    Vector3D point3 = Vector3D::point(-1, 0, 0);
    Vector3D point4 = Vector3D::point(0, -1, 0);
    Vector3D point5 = Vector3D::point(0, 0, -1);
    Vector3D point6 = Vector3D::point(0, 0, 1);
    std::vector<Vector3D> octahedron_points = {point1, point2, point3, point4, point5, point6};
    Face line1({0, 1, 5});
    Face line2({1, 2, 5});
    Face line3({2, 3, 5});
    Face line4({3, 0, 5});
    Face line5({1, 0, 4});
    Face line6({2, 1, 4});
    Face line7({3, 2, 4});
    Face line8({0, 3, 4});
    std::vector<Face> octahedron_lines = {line1, line2, line3, line4, line5, line6, line7, line8};
    figure3D figure;
    figure.points = octahedron_points;
    figure.faces = octahedron_lines;
    return figure;
}

figure3D figure3D::make_icosahedron() {
    Vector3D point1 = Vector3D::point(0, 0, sqrt(5)/2);
    std::vector<Vector3D> icosahedron_points = {point1};
    for (int i = 0; i < 5; ++i) {
        Vector3D point2_3_4_5_6 = Vector3D::point(cos(((i-2)*2*M_PI)/5), sin(((i-2)*2*M_PI)/5), 0.5);
        icosahedron_points.push_back(point2_3_4_5_6);
    }
    for (int i = 0; i < 5; ++i) {
        Vector3D point7_8_9_10_11 = Vector3D::point(cos(M_PI/5 + ((i-7)*2*M_PI)/5), sin(M_PI/5 + ((i-7)*2*M_PI)/5), -0.5);
        icosahedron_points.push_back(point7_8_9_10_11);
    }
    Vector3D point12 = Vector3D::point(0, 0, -sqrt(5)/2);
    icosahedron_points.push_back(point12);
    Face line1({0, 1, 2});
    Face line2({0, 2, 3});
    Face line3({0, 3, 4});
    Face line4({0, 4, 5});
    Face line5({0, 5, 1});
    Face line6({1, 6, 2});
    Face line7({2, 6, 7});
    Face line8({2, 7, 3});
    Face line9({3, 7, 8});
    Face line10({3, 8, 4});
    Face line11({4, 8, 9});
    Face line12({4, 9, 5});
    Face line13({5, 9, 10});
    Face line14({5, 10, 1});
    Face line15({1, 10, 6});
    Face line16({11, 7, 6});
    Face line17({11, 8, 7});
    Face line18({11, 9, 8});
    Face line19({11, 10, 9});
    Face line20({11, 6, 10});
    std::vector<Face> icosahedron_lines = {line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11,
                                           line12, line13, line14, line15, line16, line17, line18, line19, line20};
    figure3D figure;
    figure.points = icosahedron_points;
    figure.faces = icosahedron_lines;
    return figure;
}

figure3D figure3D::make_dodecahedron() {
    auto help_icosahedron = make_icosahedron();
    std::vector<Vector3D> dodecahedron_points;
    for (auto & face : help_icosahedron.faces) {
        Vector3D dodeca_point;
        dodeca_point.x = (help_icosahedron.points[face.point_indexes[0]].x + help_icosahedron.points[face.point_indexes[1]].x + help_icosahedron.points[face.point_indexes[2]].x)/3;
        dodeca_point.y = (help_icosahedron.points[face.point_indexes[0]].y + help_icosahedron.points[face.point_indexes[1]].y + help_icosahedron.points[face.point_indexes[2]].y)/3;
        dodeca_point.z = (help_icosahedron.points[face.point_indexes[0]].z + help_icosahedron.points[face.point_indexes[1]].z + help_icosahedron.points[face.point_indexes[2]].z)/3;
        dodecahedron_points.push_back(dodeca_point);
    }
    Face line1({0, 1, 2, 3, 4});
    Face line2({0, 5, 6, 7, 1});
    Face line3({1, 7, 8, 9, 2});
    Face line4({2, 9, 10, 11, 3});
    Face line5({3, 11, 12, 13, 4});
    Face line6({4, 13, 14, 5, 0});
    Face line7({19, 18, 17, 16, 15});
    Face line8({19, 14, 13, 12, 18});
    Face line9({18, 12, 11, 10, 17});
    Face line10({17, 10, 9, 8, 16});
    Face line11({16, 8, 7, 6, 15});
    Face line12({15, 6, 5, 14, 19});
    std::vector<Face> dodecahedron_lines = {line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12};
    figure3D figure;
    figure.points = dodecahedron_points;
    figure.faces = dodecahedron_lines;
    return figure;
}

figure3D figure3D::make_buckyball() {
    //TODO:: Nog niet af
    auto help_icosahedron = make_icosahedron();
    figure3D bucky;
    for(auto face: help_icosahedron.faces) {
        for(int i=0; i < face.point_indexes.size()-1; ++i) {
            //bucky.points.push_back(help_icosahedron.points[face.point_indexes[i]]);
            Vector3D firstnew_point = help_icosahedron.points[face.point_indexes[i]] + (help_icosahedron.points[face.point_indexes[i+1]] - help_icosahedron.points[face.point_indexes[i]]) * 1/3;
            bucky.points.push_back(firstnew_point);
            Vector3D scndnew_point = help_icosahedron.points[face.point_indexes[i]] + (help_icosahedron.points[face.point_indexes[i+1]] - help_icosahedron.points[face.point_indexes[i]]) * 2/3;
            bucky.points.push_back(scndnew_point);
        }
        //bucky.points.push_back(help_icosahedron.points[face.point_indexes.back()]);
        Vector3D firstnew_point = help_icosahedron.points[face.point_indexes.back()] + (help_icosahedron.points[face.point_indexes[0]] - help_icosahedron.points[face.point_indexes.back()]) * 1/3;
        bucky.points.push_back(firstnew_point);
        Vector3D scndnew_point = help_icosahedron.points[face.point_indexes.back()] + (help_icosahedron.points[face.point_indexes[0]] - help_icosahedron.points[face.point_indexes.back()]) * 2/3;
        bucky.points.push_back(scndnew_point);

        Face new_face;
        new_face.point_indexes.push_back(bucky.points.size()-6);
        new_face.point_indexes.push_back(bucky.points.size()-5);
        new_face.point_indexes.push_back(bucky.points.size()-4);
        new_face.point_indexes.push_back(bucky.points.size()-3);
        new_face.point_indexes.push_back(bucky.points.size()-2);
        new_face.point_indexes.push_back(bucky.points.size()-1);
        bucky.faces.push_back(new_face);
    }
    for(auto& point_1: help_icosahedron.points) {
        std::vector<double> vijfhoek;
        for(auto& point_2: bucky.points) {
            Vector3D line = point_2 - point_1;
            double length = line.length();
            vijfhoek.emplace_back(length);
        }
        std::vector<int> hulp;
        for(int i = 0; i < 5; ++i) {
            double val = *min_element(vijfhoek.begin(), vijfhoek.end());
            auto it = std::find(vijfhoek.begin(), vijfhoek.end(), val);
            hulp.push_back(std::distance(vijfhoek.begin(), it));
            vijfhoek.erase(vijfhoek.begin() + std::distance(vijfhoek.begin(), it));
        }
        Face new_face;
        new_face.point_indexes = hulp;
        bucky.faces.push_back(new_face);
    }

    /*
    Face new_face;
    new_face.point_indexes = {0,5,10,15,20};
    bucky.faces.push_back(new_face);
    */
    return bucky;
}

figure3D figure3D::make_cone( const int n, const double height) {
    Vector3D point1 = Vector3D::point(0, 0, height);
    std::vector<Vector3D> cone_points = {point1};
    Face linen;
    for (int i = 0; i < n; ++i) {
        Vector3D pointn = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), 0);
        cone_points.push_back(pointn);
    }
    for (int i = cone_points.size()-1; i > 0; --i) {
        linen.point_indexes.push_back(i);
    }
    std::vector<Face> cone_lines = {linen};
    for (int i = 1; i < cone_points.size()-1; ++i) {
        Face mantle;
        mantle.point_indexes.push_back(i);
        mantle.point_indexes.push_back(i+1);
        mantle.point_indexes.push_back(0);
        cone_lines.push_back(mantle);
    }
    Face mantle;
    mantle.point_indexes.push_back(cone_points.size()-1);
    mantle.point_indexes.push_back(1);
    mantle.point_indexes.push_back(0);
    cone_lines.push_back(mantle);
    figure3D figure;
    figure.points = cone_points;
    figure.faces = cone_lines;
    return figure;
}

figure3D figure3D::make_cylinder(const int n, const double height) {
    std::vector<Vector3D> cylinder_points;
    Face linen_ground;
    for (int i = 0; i < n; ++i) {
        Vector3D pointn = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), 0);
        cylinder_points.push_back(pointn);
    }
    for(int i = 0; i < cylinder_points.size(); ++i) {
        linen_ground.point_indexes.push_back(i);
    }
    std::vector<Face> cylinder_lines = {linen_ground};
    Face linen_roof;
    for (int i = 0; i < n; ++i) {
        Vector3D pointn = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), height);
        cylinder_points.push_back(pointn);
    }
    for(int i = n; i < cylinder_points.size(); ++i) {
        linen_roof.point_indexes.push_back(i);
    }
    cylinder_lines.push_back(linen_roof);
    for (int i = 0; i < n-1; ++i) {
        Face mantle;
        mantle.point_indexes.push_back(i);
        mantle.point_indexes.push_back(i+1);
        mantle.point_indexes.push_back(n+(i+1));
        mantle.point_indexes.push_back(n+i);
        cylinder_lines.push_back(mantle);
    }
    Face mantle;
    mantle.point_indexes.push_back(n-1);
    mantle.point_indexes.push_back(0);
    mantle.point_indexes.push_back(n);
    mantle.point_indexes.push_back(2*n-1);
    cylinder_lines.push_back(mantle);
    figure3D figure;
    figure.points = cylinder_points;
    figure.faces = cylinder_lines;
    return figure;
}

figure3D figure3D::make_sphere(int n) {
    // 1) Genereer een icosahedron
    auto help_icosahedron = make_icosahedron();
    std::vector<Vector3D> sphere_points;
    std::vector<Face> sphere_lines;
    // 2) Deel elke driehoek op  in kleinere driehoeken
    // 3) Herhaal stap 2 in totaal ‘n’ keer
    while(n != 0) {
        for (auto &face : help_icosahedron.faces) {
            //Face line;
            for (int i = 0; i < face.point_indexes.size()-1; ++i) {
                sphere_points.push_back(Vector3D::normalise(help_icosahedron.points[face.point_indexes[i]]));
                auto a = (help_icosahedron.points[face.point_indexes[i]] + help_icosahedron.points[face.point_indexes[i+1]]) / 2;
                a.normalise();
                sphere_points.push_back(a);
            }
            sphere_points.push_back(Vector3D::normalise(help_icosahedron.points[face.point_indexes[face.point_indexes.size()-1]]));
            auto a = (help_icosahedron.points[face.point_indexes[face.point_indexes.size()-1]] + help_icosahedron.points[face.point_indexes[0]]) / 2;
            a.normalise();
            sphere_points.push_back(a);

            int place = sphere_points.size();
            Face rectangle1({place-6, place-5, place-1});
            Face rectangle2({place-5, place-4, place-3});
            Face rectangle3({place-3, place-2, place-1});
            Face rectangle4({place-5, place-3, place-1});
            sphere_lines.push_back(rectangle1);
            sphere_lines.push_back(rectangle2);
            sphere_lines.push_back(rectangle3);
            sphere_lines.push_back(rectangle4);
        }
        figure3D figure;
        figure.points = sphere_points;
        figure.faces = sphere_lines;
        sphere_lines.clear();
        help_icosahedron = figure;
        --n;
    }
    return help_icosahedron;
}

figure3D figure3D::make_torus(const double r, const double R, const int n, const int m) {
    std::vector<Vector3D> torus_points;
    for (int i = n-1; i >= 0; --i) {
        for (int j = m-1; j >= 0; --j) {
            double u = (2*i*M_PI)/n;
            double v = (2*j*M_PI)/m;
            Vector3D point = Vector3D::point((R+(r * cos(v))) * cos(u), (R+(r * cos(v))) * sin(u), r*sin(v));
            torus_points.push_back(point);
        }
    }
    std::vector<Face> torus_lines;
    for (int i = n-1; i >= 0; --i) {
        for (int j = m-1; j >= 0; --j) {
            Face line;
            line.point_indexes = {i*m +j, ((i+1)%n)*m+j, ((i+1)%n)*m+((j+1)%m), i*m + ((j+1)%m)};
            torus_lines.push_back(line);
        }
    }
    figure3D figure;
    figure.points = torus_points;
    figure.faces = torus_lines;
    return figure;
}

figure3D translate_cube(figure3D fig, Matrix frac_scale, int current_it, double x, double y, double z) {
    figure3D new_figure = fig;
    for(auto& smaller : new_figure.points) {
        smaller *= frac_scale;
    }
    Matrix T;
    T(4,1) = x;
    T(4,2) = y;
    T(4,3) = z;
    for(auto& smaller : new_figure.points) {
        smaller *= T;
    }
    return new_figure;
}

std::list<figure3D> figure3D::make_menger(figure3D& fig, std::list<figure3D> fractals, const int nr_iterations, int current_it) {
    if(current_it == nr_iterations) {
        return fractals;
    }
    Matrix frac_scale;
    Figures3D temps;
    frac_scale(1, 1) = 1/(std::pow(3, current_it));
    frac_scale(2, 2) = 1/(std::pow(3, current_it));
    frac_scale(3, 3) = 1/(std::pow(3, current_it));

    for(auto figure : fractals) {
        for (int i = 0; i < figure.points.size(); ++i) {
            figure3D new_figure = fig;
            for(auto& smaller : new_figure.points) {
                smaller *= frac_scale;
            }
            Matrix T;
            T(4,1) = figure.points[i].x - new_figure.points[i].x;
            T(4,2) = figure.points[i].y - new_figure.points[i].y;
            T(4,3) = figure.points[i].z - new_figure.points[i].z;
            for(auto& smaller : new_figure.points) {
                smaller *= T;
            }
            temps.push_back(new_figure);
        }

        // Dit is voor laag 1
        temps.push_back(translate_cube(fig, frac_scale, current_it, 2/(std::pow(3, current_it)), 0, -2/(std::pow(3, current_it))));
        temps.push_back(translate_cube(fig, frac_scale, current_it, -2/(std::pow(3, current_it)), 0, -2/(std::pow(3, current_it))));
        temps.push_back(translate_cube(fig, frac_scale, current_it, 0, 2/(std::pow(3, current_it)), -2/(std::pow(3, current_it))));
        temps.push_back(translate_cube(fig, frac_scale, current_it, 0, -2/(std::pow(3, current_it)), -2/(std::pow(3, current_it))));

        // Dit is voor diegene waar de z-as niet voor verandert, laag 2
        temps.push_back(translate_cube(fig, frac_scale, current_it, 2/(std::pow(3, current_it)), -2/(std::pow(3, current_it)), 0));
        temps.push_back(translate_cube(fig, frac_scale, current_it, -2/(std::pow(3, current_it)), -2/(std::pow(3, current_it)), 0));
        temps.push_back(translate_cube(fig, frac_scale, current_it, 2/(std::pow(3, current_it)), 2/(std::pow(3, current_it)), 0));
        temps.push_back(translate_cube(fig, frac_scale, current_it, -2/(std::pow(3, current_it)), 2/(std::pow(3, current_it)), 0));

        // Dit is voor laag 3
        temps.push_back(translate_cube(fig, frac_scale, current_it, 2/(std::pow(3, current_it)), 0, 2/(std::pow(3, current_it))));
        temps.push_back(translate_cube(fig, frac_scale, current_it, -2/(std::pow(3, current_it)), 0, 2/(std::pow(3, current_it))));
        temps.push_back(translate_cube(fig, frac_scale, current_it, 0, 2/(std::pow(3, current_it)), 2/(std::pow(3, current_it))));
        temps.push_back(translate_cube(fig, frac_scale, current_it, 0, -2/(std::pow(3, current_it)), 2/(std::pow(3, current_it))));

    }
    fractals = temps;
    return make_menger(fig, fractals, nr_iterations, current_it+=1);
}

std::string replaceD3L(const std::string& replacement, int iterations, const LParser::LSystem3D& l_system, std::string total){
    for (char i : replacement) {
        if(i == '+' or i == '-' or i == '^' or i == '&' or i == '\\' or i == '/' or i == '(' or i == ')') {
            total.push_back(i);
        }
        else{
            total += (l_system.get_replacement(i));
        }
    }
    while (iterations != 1) {
        iterations -= 1;
        return replaceD3L(total, iterations, l_system, "");
    }
    return total;
}

figure3D figure3D::draw3LSystem(std::string inputfile) {
    // Zet alles klaar voor de l_parser
    LParser::LSystem3D l_system;
    std::ifstream input_stream(inputfile);
    input_stream >> l_system;
    input_stream.close();

    // Haal alle informatie uit de.L3D file
    auto alfabet = l_system.get_alphabet();
    unsigned int iterations = l_system.get_nr_iterations();
    std::string initiator = l_system.get_initiator();
    double addangle = l_system.get_angle();
    std::string bigstring;

    // Grote string maken voor de lijnen te tekenen
    for (char i : initiator) {
        if(i == '+' or i == '-' or i == '^' or i == '&' or i == '\\' or i == '/') {
            bigstring += i;
        } else {
            bigstring += replaceD3L(l_system.get_replacement(i), iterations-1, l_system, "");
        }
    }
    // Grote string lezen en figuur maken
    figure3D L3D;
    Vector3D currentpoint;
    std::vector<std::tuple<Vector3D, Vector3D, Vector3D, Vector3D>> stack;

    Vector3D H_new = Vector3D::vector(1,0,0);
    Vector3D L_new = Vector3D::vector(0,1,0);
    Vector3D U_new = Vector3D::vector(0,0,1);
    Vector3D H = Vector3D::vector(1,0,0);
    Vector3D L = Vector3D::vector(0,1,0);
    Vector3D U = Vector3D::vector(0,0,1);

    for (char j : bigstring) {
        if (j == '+') {
            H_new = H * cos(addangle * (M_PI/180)) + L * sin(addangle * (M_PI/180));
            L_new = -H * sin(addangle * (M_PI/180)) + L * cos(addangle * (M_PI/180));
        }
        else if (j == '-') {
            H_new = H * cos(-addangle * (M_PI/180)) + L * sin(-addangle * (M_PI/180));
            L_new = -H * sin(-addangle * (M_PI/180)) + L * cos(-addangle * (M_PI/180));
        }
        else if (j == '(') {
            stack.emplace_back(currentpoint, H, L, U);
        }
        else if (j == ')') {
            currentpoint = std::get<0>(stack.back());
            H_new = std::get<1>(stack.back());
            L_new = std::get<2>(stack.back());
            U_new = std::get<3>(stack.back());
            stack.pop_back();
        }
        else if (j == '^') {
            H_new = H * cos(addangle * (M_PI/180)) + U * sin(addangle * (M_PI/180));
            U_new = -H * sin(addangle * (M_PI/180)) + U * cos(addangle * (M_PI/180));
        }
        else if (j == '&') {
            H_new = H * cos(-addangle * (M_PI/180)) + U * sin(-addangle * (M_PI/180));
            U_new = -H * sin(-addangle * (M_PI/180)) + U * cos(-addangle * (M_PI/180));
        }
        else if (j == '\\') {
            L_new = L * cos(addangle * (M_PI/180)) - U * sin(addangle * (M_PI/180));
            U_new = L * sin(addangle * (M_PI/180)) + U * cos(addangle * (M_PI/180));
        }
        else if (j == '/') {
            L_new = L * cos(-addangle * (M_PI/180)) - U * sin(-addangle * (M_PI/180));
            U_new = L * sin(-addangle * (M_PI/180)) + U * cos(-addangle * (M_PI/180));
        }
        else {
            Vector3D new_point;
            new_point = currentpoint + H;
            if(l_system.draw(j)){
                L3D.points.push_back(currentpoint);
                L3D.points.push_back(new_point);
                Face line;
                line.point_indexes.push_back(L3D.points.size()-2);
                line.point_indexes.push_back(L3D.points.size()-1);
                L3D.faces.push_back(line);
            }
            currentpoint = new_point;
        }
        H = H_new;
        L = L_new;
        U = U_new;
    }
    return L3D;
}

Matrix scaleFigure(const double scale) {
    Matrix scaler;
    for (int i = 1; i < 4; ++i) {
        scaler(i, i) = scale;
    }
    return scaler;
}

Matrix rotateX(double angle) {
    angle *= M_PI/180;
    Matrix rotator;
    rotator(2,2) = cos(angle);
    rotator(2, 3) = sin(angle);
    rotator(3,3) = cos(angle);
    rotator(3, 2) = -sin(angle);
    return rotator;
}

Matrix rotateY(double angle) {
    angle *= M_PI/180;
    Matrix rotator;
    rotator(1,1) = cos(angle);
    rotator(1, 3) = -sin(angle);
    rotator(3,1) = sin(angle);
    rotator(3, 3) = cos(angle);
    return rotator;
}

Matrix rotateZ(double angle) {
    angle *= M_PI/180;
    Matrix rotator;
    rotator(1,1) = cos(angle);
    rotator(1,2) = sin(angle);
    rotator(2, 1) = -sin(angle);
    rotator(2, 2) = cos(angle);
    return rotator;
}

Matrix translate(const Vector3D &vector) {
    Matrix translator;
    translator(4, 1) = vector.x;
    translator(4, 2) = vector.y;
    translator(4, 3) = vector.z;
    return translator;
}

void applyTransformation(figure3D &figure, const Matrix &matrix) {
    for (auto & point : figure.points) {
        point *= matrix;
    }
}

Matrix eyePointTrans(const Vector3D &eyepoint) {
    Matrix eye;
    double theta = 0;
    double phi = 0;
    double r = 0;
    toPolar(eyepoint, theta, phi, r);
    eye(1,1) = -sin(theta);
    eye(1,2) = -cos(theta) * cos(phi);
    eye(1,3) = cos(theta) * sin(phi);
    eye(2,1) = cos(theta);
    eye(2,2) = -sin(theta) * cos(phi);
    eye(2,3) = sin(theta) * sin(phi);
    eye(3, 2) = sin(phi);
    eye(3,3) = cos(phi);
    eye(4, 3) = -r;
    return eye;
}

void toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    // Omzetten van carthesische naar poolcoördinaten
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = atan2(point.y, point.x);
    phi = acos(point.z / r);
}

Lines2D doProjection(const Figures3D &figures) {
    std::list<Line2D> lines;
    for(auto i: figures) {
        std::vector<Point2D> points;
        for(auto j: i.points) {
            points.push_back(doProjection(j, 1));
        }
        for(auto j: i.faces) {
            Line2D line(points[j.point_indexes.back()], points[j.point_indexes[0]], i.ambientReflection, i.points[j.point_indexes.back()].z, i.points[j.point_indexes[0]].z);
            lines.push_back(line);
            for (int k = 0; k < j.point_indexes.size()-1; ++k) {
                Line2D line(points[j.point_indexes[k]], points[j.point_indexes[k+1]], i.ambientReflection, i.points[j.point_indexes[k]].z, i.points[j.point_indexes[k+1]].z);
                lines.push_back(line);
            }
        }
    }
    return lines;
}

void applyTransformation(Figures3D &figures, const Matrix &matrix) {
    for (auto& i: figures) {
        applyTransformation(i, matrix);
    }
}

Point2D doProjection(const Vector3D &point, const double d) {
    Point2D newpoint;
    newpoint.x = -(point.x*d / point.z);
    newpoint.y = -(point.y*d / point.z);
    return newpoint;
}

img::EasyImage line_drawings(const ini::Configuration &configuration) {
    int number_of_figures = configuration["General"]["nrFigures"].as_int_or_die();
    std::vector<double> eyevec = configuration["General"]["eye"].as_double_tuple_or_die();
    int size = configuration["General"]["size"].as_int_or_die();
    std::vector<double> background = configuration["General"]["backgroundcolor"];
    Color backgr(background[0], background[1], background[2]);
    img::Color back = convertcolor(backgr);
    Vector3D eye = Vector3D::point(eyevec[0], eyevec[1], eyevec[2]);
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    std::list<figure3D> figures;
    for (int i = 0; i < number_of_figures; ++i) {
        std::vector<double> foreground;
        bool ambient = configuration["Figure" + std::to_string(i)]["ambientReflection"].as_double_tuple_if_exists(foreground);
        if(!ambient) {
            foreground = configuration["Figure" + std::to_string(i)]["color"];
        }
        Color fore(foreground[0], foreground[1], foreground[2]);

        std::vector<double> diffusevector;
        bool diffuse = configuration["Figure" + std::to_string(i)]["diffuseReflection"].as_double_tuple_if_exists(diffusevector);
        Color diffusecolor = Color(0, 0, 0);
        if(diffuse) {
            diffusecolor = Color(diffusevector[0], diffusevector[1], diffusevector[2]);
        }

        std::vector<double> specularvector;
        bool specular = configuration["Figure" + std::to_string(i)]["specularReflection"].as_double_tuple_if_exists(specularvector);
        Color specularcolor = Color(0, 0, 0);
        if(specular) {
            specularcolor = Color(specularvector[0], specularvector[1], specularvector[2]);
        }

        double refl_coefficient = 0;
        configuration["Figure" + std::to_string(i)]["reflectionCoefficient"].as_double_if_exists(refl_coefficient);

        std::vector<double> centre = configuration["Figure" + std::to_string(i)]["center"].as_double_tuple_or_die();
        Vector3D center;
        center.x = centre[0];
        center.y = centre[1];
        center.z = centre[2];
        double angleX = configuration["Figure" + std::to_string(i)]["rotateX"].as_double_or_die();
        double angleY = configuration["Figure" + std::to_string(i)]["rotateY"].as_double_or_die();
        double angleZ = configuration["Figure" + std::to_string(i)]["rotateZ"].as_double_or_die();
        double scale = configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die();
        Matrix M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
        figure3D figure;
        figure.ambientReflection = fore;
        figure.diffuseReflection = diffusecolor;
        figure.specularReflection = specularcolor;
        figure.reflectionCoefficient = refl_coefficient;

        if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "LineDrawing") {
            int nrPoints = configuration["Figure" + std::to_string(i)]["nrPoints"].as_int_or_die();
            int nrLines = configuration["Figure" + std::to_string(i)]["nrLines"].as_int_or_die();
            for (int j = 0; j < nrPoints; ++j) {
                //Vector3D point
                std::vector<double> info = configuration["Figure" + std::to_string(i)]["point" + std::to_string(j)];
                Vector3D point = Vector3D::point(info[0], info[1], info[2]);
                points.push_back(point);
            }
            for (int j = 0; j < nrLines; ++j) {
                //Face line
                std::vector<double> info = configuration["Figure" + std::to_string(i)]["line" + std::to_string(j)];
                Face line;
                line.point_indexes.push_back(info[0]);
                line.point_indexes.push_back(info[1]);
                faces.push_back(line);
            }
            figure.points = points;
            figure.faces = faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cube") {
            figure3D temp_fig = figure3D::make_cube();
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Tetrahedron") {
            figure3D temp_fig = figure3D::make_tetrahedron();
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Octahedron") {
            figure3D temp_fig = figure3D::make_octahedron();
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Icosahedron") {
            figure3D temp_fig = figure3D::make_icosahedron();
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Dodecahedron") {
            figure3D temp_fig = figure3D::make_dodecahedron();
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Sphere") {
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            figure3D temp_fig = figure3D::make_sphere(n);
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cone") {
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
            figure3D temp_fig = figure3D::make_cone(n, height);
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cylinder") {
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
            figure3D temp_fig = figure3D::make_cylinder(n, height);
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Torus") {
            double r = configuration["Figure" + std::to_string(i)]["r"].as_double_or_die();
            double R = configuration["Figure" + std::to_string(i)]["R"].as_double_or_die();
            int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            int m = configuration["Figure" + std::to_string(i)]["m"].as_int_or_die();
            figure3D temp_fig = figure3D::make_torus(r, R, n, m);;
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "BuckyBall") {
            figure3D temp_fig = figure3D::make_buckyball();
            figure.points = temp_fig.points;
            figure.faces = temp_fig.faces;
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "MengerSponge") {
            //TODO ;; als af kijk nog voor licht aan teoassen hier
            int iterations = configuration["Figure" + std::to_string(i)]["nrIterations"].as_int_or_die();
            figure3D hulp_cube = figure3D::make_cube();
            Figures3D fractals = {figure};
            fractals = figure3D::make_menger(hulp_cube, {hulp_cube}, iterations+1, 1);
            for(auto fractal : fractals) {
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(fractal, M);
                M = eyePointTrans(eye);
                applyTransformation(fractal, M);
                figures.push_back(fractal);
            }
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "3DLSystem") {
            std::string inputfile = configuration["Figure" + std::to_string(i)]["inputfile"].as_string_or_die();
            figure = figure3D::draw3LSystem(inputfile);
            applyTransformation(figure, M);
            M = eyePointTrans(eye);
            applyTransformation(figure, M);
            figures.push_back(figure);
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die().find("Fractal") != std::string::npos) {
            int iterations = configuration["Figure" + std::to_string(i)]["nrIterations"].as_int_or_die();
            double frac_scale = configuration["Figure" + std::to_string(i)]["fractalScale"].as_double_or_die();

            if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalTetrahedron") {
                figure3D temp_fig = figure3D::make_tetrahedron();
                figure.points = temp_fig.points;
                figure.faces = temp_fig.faces;
            }
            else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalCube") {
                figure3D temp_fig = figure3D::make_cube();
                figure.points = temp_fig.points;
                figure.faces = temp_fig.faces;
            }
            else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalIcosahedron") {
                figure3D temp_fig = figure3D::make_icosahedron();
                figure.points = temp_fig.points;
                figure.faces = temp_fig.faces;
            }
            else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalOctahedron") {
                figure3D temp_fig = figure3D::make_octahedron();
                figure.points = temp_fig.points;
                figure.faces = temp_fig.faces;
            }
            else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalDodecahedron") {
                figure3D temp_fig = figure3D::make_dodecahedron();
                figure.points = temp_fig.points;
                figure.faces = temp_fig.faces;
            }
            Figures3D fractals = {figure};
            fractals = generateFractal(figure, fractals, iterations+1, 1,  frac_scale);
            for(auto fractal : fractals) {
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(fractal, M);
                M = eyePointTrans(eye);
                applyTransformation(fractal, M);
                figures.push_back(fractal);
            }
        }

        // THICKKKKK
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() ==  "ThickLineDrawing") {
            int nrPoints = configuration["Figure" + std::to_string(i)]["nrPoints"].as_int_or_die();
            int nrLines = configuration["Figure" + std::to_string(i)]["nrLines"].as_int_or_die();
            double r = configuration["Figure" + std::to_string(i)]["radius"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            int m = configuration["Figure" + std::to_string(i)]["m"];
            for (int j = 0; j < nrPoints; ++j) {
                std::vector<double> info = configuration["Figure" + std::to_string(i)]["point" + std::to_string(j)];
                Vector3D point = Vector3D::point(info[0], info[1], info[2]);
                points.push_back(point);
            }
            for (int j = 0; j < nrLines; ++j) {
                std::vector<double> info = configuration["Figure" + std::to_string(i)]["line" + std::to_string(j)];
                Face line;
                line.point_indexes.push_back(info[0]);
                line.point_indexes.push_back(info[1]);
                faces.push_back(line);
            }
            figure.points = points;
            figure.faces = faces;
            Figures3D thicks = {};
            generateThickFigure(figure, thicks, r, n, m);
            for(auto thick : thicks) {
                thick.ambientReflection = figure.ambientReflection;
                thick.diffuseReflection = figure.diffuseReflection;
                thick.specularReflection = figure.specularReflection;
                thick.reflectionCoefficient = figure.reflectionCoefficient;
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(thick, M);
                M = eyePointTrans(eye);
                applyTransformation(thick, M);
                figures.push_back(thick);
            }
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Thick3DLSystem") {
            std::string inputfile = configuration["Figure" + std::to_string(i)]["inputfile"].as_string_or_die();
            double r = configuration["Figure" + std::to_string(i)]["radius"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            int m = configuration["Figure" + std::to_string(i)]["m"];
            figure3D temp_figure = figure3D::draw3LSystem(inputfile);
            figure.points = temp_figure.points;
            figure.faces = temp_figure.faces;
            Figures3D thicks = {};
            generateThickFigure(figure, thicks, r, n, m);
            for(auto thick : thicks) {
                thick.ambientReflection = figure.ambientReflection;
                thick.diffuseReflection = figure.diffuseReflection;
                thick.specularReflection = figure.specularReflection;
                thick.reflectionCoefficient = figure.reflectionCoefficient;
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(thick, M);
                M = eyePointTrans(eye);
                applyTransformation(thick, M);
                figures.push_back(thick);
            }
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "ThickCube") {
            double r = configuration["Figure" + std::to_string(i)]["radius"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            int m = configuration["Figure" + std::to_string(i)]["m"];
            figure3D temp_fig = figure3D::make_cube();
            Figures3D thicks = {};
            generateThickFigure(temp_fig, thicks, r, n, m);
            for(auto thick : thicks) {
                thick.ambientReflection = figure.ambientReflection;
                thick.diffuseReflection = figure.diffuseReflection;
                thick.specularReflection = figure.specularReflection;
                thick.reflectionCoefficient = figure.reflectionCoefficient;
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(thick, M);
                M = eyePointTrans(eye);
                applyTransformation(thick, M);
                figures.push_back(thick);
            }
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "ThickTetrahedron") {
            double r = configuration["Figure" + std::to_string(i)]["radius"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            int m = configuration["Figure" + std::to_string(i)]["m"];
            figure3D temp_fig = figure3D::make_tetrahedron();
            Figures3D thicks = {};
            generateThickFigure(temp_fig, thicks, r, n, m);
            for(auto thick : thicks) {
                thick.ambientReflection = figure.ambientReflection;
                thick.diffuseReflection = figure.diffuseReflection;
                thick.specularReflection = figure.specularReflection;
                thick.reflectionCoefficient = figure.reflectionCoefficient;
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(thick, M);
                M = eyePointTrans(eye);
                applyTransformation(thick, M);
                figures.push_back(thick);
            }
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "ThickIcosahedron") {
            double r = configuration["Figure" + std::to_string(i)]["radius"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            int m = configuration["Figure" + std::to_string(i)]["m"];
            figure3D temp_fig = figure3D::make_icosahedron();
            Figures3D thicks = {};
            generateThickFigure(temp_fig, thicks, r, n, m);
            for(auto thick : thicks) {
                thick.ambientReflection = figure.ambientReflection;
                thick.diffuseReflection = figure.diffuseReflection;
                thick.specularReflection = figure.specularReflection;
                thick.reflectionCoefficient = figure.reflectionCoefficient;
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(thick, M);
                M = eyePointTrans(eye);
                applyTransformation(thick, M);
                figures.push_back(thick);
            }
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "ThickOctahedron") {
            double r = configuration["Figure" + std::to_string(i)]["radius"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            int m = configuration["Figure" + std::to_string(i)]["m"];
            figure3D temp_fig = figure3D::make_octahedron();
            Figures3D thicks = {};
            generateThickFigure(temp_fig, thicks, r, n, m);
            for(auto thick : thicks) {
                thick.ambientReflection = figure.ambientReflection;
                thick.diffuseReflection = figure.diffuseReflection;
                thick.specularReflection = figure.specularReflection;
                thick.reflectionCoefficient = figure.reflectionCoefficient;
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(thick, M);
                M = eyePointTrans(eye);
                applyTransformation(thick, M);
                figures.push_back(thick);
            }
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "ThickDodecahedron") {
            double r = configuration["Figure" + std::to_string(i)]["radius"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            int m = configuration["Figure" + std::to_string(i)]["m"];
            figure3D temp_fig = figure3D::make_dodecahedron();
            Figures3D thicks = {};
            generateThickFigure(temp_fig, thicks, r, n, m);
            for(auto thick : thicks) {
                thick.ambientReflection = figure.ambientReflection;
                thick.diffuseReflection = figure.diffuseReflection;
                thick.specularReflection = figure.specularReflection;
                thick.reflectionCoefficient = figure.reflectionCoefficient;
                M = scaleFigure(scale) * rotateX(angleX) * rotateY(angleY) * rotateZ(angleZ) * translate(center);
                applyTransformation(thick, M);
                M = eyePointTrans(eye);
                applyTransformation(thick, M);
                figures.push_back(thick);
            }
        }
    }

    int nrLights = configuration["General"]["nrLights"].as_int_or_default(0);
    Lights3D lightings;
    for (int i = 0; i < nrLights; ++i) {
        std::vector<double> ambients = configuration["Light" + std::to_string(i)]["ambientLight"];
        bool infinity;
        if(configuration["Light" + std::to_string(i)]["infinity"].as_bool_if_exists(infinity)) {
            std::vector<double> diffuses = configuration["Light" + std::to_string(i)]["diffuseLight"];
            std::vector<double> speculars = configuration["Light" + std::to_string(i)]["specularLight"].as_double_tuple_or_default({0, 0, 0});
            if(infinity) {
                std::vector<double> direction = configuration["Light" + std::to_string(i)]["direction"].as_double_tuple_or_die();
                auto* inflight = new InfLight(Color(ambients[0], ambients[1], ambients[2]), Color(diffuses[0], diffuses[1], diffuses[2]), Color(speculars[0],speculars[1],speculars[2]), Vector3D::vector(direction[0], direction[1], direction[2]));
                inflight->ldVector *= eyePointTrans(eye);
                inflight->eye_trans = eye*eyePointTrans(eye);
                lightings.push_back(inflight);
            } else {
                std::vector<double> location = configuration["Light" + std::to_string(i)]["location"].as_double_tuple_or_die();
                Vector3D locationpoint = Vector3D::point(location[0], location[1], location[2]);
                locationpoint *= eyePointTrans(eye);
                double spotangle = configuration["Light" + std::to_string(i)]["spotAngle"].as_double_or_default(90);
                auto* pointlight = new PointLight(Color(ambients[0], ambients[1], ambients[2]),
                                                  Color(diffuses[0], diffuses[1], diffuses[2]),
                                                  Color(speculars[0], speculars[1], speculars[2]), locationpoint,
                                                  spotangle);
                pointlight->eye_trans = eye*eyePointTrans(eye);
                lightings.push_back(pointlight);
            }
        } else {
            auto* light = new Light(Color(ambients[0], ambients[1], ambients[2]), Color(0, 0, 0), Color(0, 0, 0));
            lightings.push_back(light);
        }
    }
    if(lightings.empty()) {
        auto* nieuw = new Light(Color(1,1,1), Color(0,0,0), Color(0,0,0));
        lightings.push_back(nieuw);
    }
    if(configuration["General"]["type"].as_string_or_die() == "ZBuffering" or configuration["General"]["type"].as_string_or_die() == "LightedZBuffering") {
        for(auto & figure : figures) {
            std::vector<Face> triangle_faces;
            for(auto & face : figure.faces) {
                auto triangles = triangulate(face);
                for(auto & triangle : triangles) {
                    triangle_faces.push_back(triangle);
                }
            }
            figure.faces = triangle_faces;
        }
        return draw3DLines(figures, size, back, lightings);
    }
    auto lines = doProjection(figures);
    return draw2DLines(lines, size, back, configuration["General"]["type"].as_string_or_die());
}

std::vector<Face> triangulate(const Face& face) {
    std::vector<Face> faces;
    for(int i = 1; i < face.point_indexes.size()-1; ++i) {
        Face triangle({face.point_indexes[0], face.point_indexes[i], face.point_indexes[i+1]});
        faces.push_back(triangle);
    }
    return faces;
}

img::EasyImage draw3DLines(Figures3D &figures, const int size, const img::Color& backgroundcolor, Lights3D lights) {
    Lines2D lines = doProjection(figures);
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

    img::EasyImage image(imagex, imagey, backgroundcolor);
    ZBuffer ZBuffer(imagex, imagey);
    for (figure3D figure:figures) {
        for(auto & face : figure.faces) {
            ZBuffer.draw_zbuf_triag(image, figure.points[face.point_indexes[0]], figure.points[face.point_indexes[1]],
            figure.points[face.point_indexes[2]], d, dx, dy, figure.ambientReflection, figure.diffuseReflection, figure.specularReflection, figure.reflectionCoefficient, lights);
        }
    }
    return image;
}

Figures3D generateFractal(figure3D& fig, Figures3D& fractals, const int nr_iterations, int current_it, double scale) {
    if(current_it == nr_iterations) {
        return fractals;
    }
    Matrix frac_scale;
    Figures3D temps;
    frac_scale(1, 1) = 1/(std::pow(scale, current_it));
    frac_scale(2, 2) = 1/(std::pow(scale, current_it));
    frac_scale(3, 3) = 1/(std::pow(scale, current_it));

    for(auto figure : fractals) {
        for (int i = 0; i < figure.points.size(); ++i) {
            figure3D new_figure = fig;
            for(auto& smaller : new_figure.points) {
                smaller *= frac_scale;
            }
            Matrix T;
            T(4,1) = figure.points[i].x - new_figure.points[i].x;
            T(4,2) = figure.points[i].y - new_figure.points[i].y;
            T(4,3) = figure.points[i].z - new_figure.points[i].z;
            for(auto& smaller : new_figure.points) {
                smaller *= T;
            }
            temps.push_back(new_figure);
        }
    }
    fractals = temps;
    return generateFractal(fig, fractals, nr_iterations, current_it+=1, scale);
}

void generateThickFigure(const figure3D &lineDrawing, Figures3D &resultingFigures, const double r, const int n, const int m) {
    // Maak van elk punt een bol
    Matrix frac_scale;
    frac_scale(1, 1) = frac_scale(2, 2) = frac_scale(3, 3) = r;
    for(auto& point: lineDrawing.points) {
        figure3D bol = figure3D::make_sphere(m);
        applyTransformation(bol, scaleFigure(r));
        applyTransformation(bol, translate(point));
        resultingFigures.push_back(bol);
    }
    // Maak van elke lijn een cillinder
    for(auto& face: lineDrawing.faces) {
        for (int i = 0; i < face.point_indexes.size(); ++i) {
            Vector3D line;
            if(i == face.point_indexes.size()-1) {
                line = lineDrawing.points[face.point_indexes[0]] - lineDrawing.points[face.point_indexes[i]];
            } else {
                line = lineDrawing.points[face.point_indexes[i+1]] - lineDrawing.points[face.point_indexes[i]];
            }
            double height = line.length()/r;
            figure3D cillinder = figure3D::make_cylinder(n, height);
            applyTransformation(cillinder, scaleFigure(r));
            double theta, phi, radius = 0;
            toPolar(line, theta, phi, radius);
            applyTransformation(cillinder, rotateY(phi * 180/M_PI) * rotateZ(theta* 180/M_PI));
            applyTransformation(cillinder, translate(lineDrawing.points[face.point_indexes[i]]));
            resultingFigures.push_back(cillinder);
        }
    }
}

figure3D::figure3D(): points(), faces(), ambientReflection(0,0,0), diffuseReflection(0,0,0), specularReflection(0, 0,0) {}

Face::Face(std::vector<int> point_indexes): point_indexes(point_indexes) {}

Face::Face() {}
