//
// Created by Sam Pieters on 06/05/2020.
//

#include "Light.h"
#include "draw.h"

Light::Light(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight) : ambientLight(
        ambientLight), diffuseLight(diffuseLight), specularLight(specularLight) {}

Light::~Light() {

}

InfLight::InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
                   const Vector3D &ldVector) : Light(ambientLight, diffuseLight, specularLight), ldVector(ldVector) {}

PointLight::PointLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
                       const Vector3D &location, double spotAngle) : Light(ambientLight, diffuseLight, specularLight),
                                                                                         location(location), spotAngle(spotAngle) {}

ZBuffer::ZBuffer(const int width, const int height) {
    for(int i = 0; i < height; ++i) {
        std::vector<double> line;
        for (int j = 0; j < width; ++j) {
            line.push_back(std::numeric_limits<double>::infinity());
        }
        this->buffer.push_back(line);
    }
}

double ZBuffer::z_vergelijking(double z_a, double z_b, int i, int a) {
    return ((double)i/(double)a)/z_b + (1-((double)i/(double)a))/z_a;
}


void ZBuffer::draw_zbuf_line(img::EasyImage & image, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, const Color &color) {
    auto a = convertcolor(color);
    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());
    if (x0 == x1) {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
            if(check_smaller(z_vergelijking(z0, z1, i, std::max(y0, y1)), i, x0)) {
                (image)(x0, i) = a;
            }
        }
    }
    else if (y0 == y1) {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            if(check_smaller(z_vergelijking(z0, z1, i, std::max(x0, x1)), y0, i)) {
                (image)(i, y0) = a;
            }
        }
    } else {
        if (x0 > x1) {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0) {
            for (unsigned int i = 0; i <= (x1-x0); i++) {
                if(check_smaller(z_vergelijking(z0, z1, i, (x1-x0)), (unsigned int) round(y0+m*i), x0+i)) {
                    (image)(x0 + i, (unsigned int) round(y0 + m * i)) = a;
                }
            }
        }
        else if (m > 1.0) {
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                if(check_smaller(z_vergelijking(z0, z1, i, (y1-y0)), y0+i, (unsigned int) round(x0+(i/m)))) {
                    (image)((unsigned int) round(x0 + (i / m)), y0 + i) = a;
                }
            }
        }
        else if (m < -1.0) {
            for (unsigned int i = 0; i <= (y0 - y1); i++) {
                if(check_smaller(z_vergelijking(z0, z1, i, (y0-y1)), y0-i, (unsigned int) round(x0-(i/m)))) {
                    (image)((unsigned int) round(x0 -(i/m)), y0-i) = a;
                }
            }
        }
    }
}


bool ZBuffer::check_smaller(double z, int row, int column) {
    if(z < this->buffer[row][column]) {
        this->buffer[row][column] = z;
        return true;
    }
    return false;
}

Color tooBright(Color total_light, Color color) {
    if(total_light.red + color.red > 1) {
        total_light.red = 1;
    } else {
        total_light.red += color.red;
    }
    if(total_light.green + color.green > 1) {
        total_light.green = 1;
    } else {
        total_light.green += color.green;
    }
    if(total_light.blue + color.blue > 1) {
        total_light.blue = 1;
    } else {
        total_light.blue += color.blue;
    }
    return total_light;
}

void ZBuffer::draw_zbuf_triag(img::EasyImage& image, Vector3D const& A, Vector3D const& B, Vector3D const& C, double d, double dx, double dy, Color ambientReflection, Color diffuseReflection, Color specularReflection, double reflectionCoeff, Lights3D& lights) {
    std::vector<Vector3D> abc = {A, B, C};
    std::vector<Point2D> proj_abc;
    Color total_light(0, 0, 0);

    // ambient light
    for(Light* light: lights) {
        Color temp = Color(ambientReflection.red * light->ambientLight.red, ambientReflection.green * light->ambientLight.green, ambientReflection.blue * light->ambientLight.blue);
        total_light.red += temp.red;
        total_light.green += temp.green;
        total_light.blue += temp.blue;
    }

    //Diffuus light
    Vector3D AB = B-A;
    Vector3D AC = C-A;
    Vector3D ABC = Vector3D::cross(AB, AC);
    Vector3D n = Vector3D::normalise(ABC);
    for(Light* light: lights) {
        if(auto* infinity = dynamic_cast<InfLight*>(light)) {
            double cosalfa = - Vector3D::dot(n, Vector3D::normalise(infinity->ldVector));
            if(cosalfa > 0) {
                Color temp = Color(diffuseReflection.red * light->diffuseLight.red * cosalfa, diffuseReflection.green * light->diffuseLight.green * cosalfa, diffuseReflection.blue * light->diffuseLight.blue * cosalfa);
                total_light = tooBright(total_light, temp);
            }
        }
    }
    for(auto & punt : abc) {
        Point2D new_point;
        new_point.x  = ((d * punt.x)/(-punt.z)) + dx;
        new_point.y = ((d * punt.y)/(-punt.z)) + dy;
        proj_abc.push_back(new_point);
    }

    int min = static_cast<int>(round(std::min(std::min(proj_abc[0].y, proj_abc[1].y), proj_abc[2].y) + 0.5));
    int max = static_cast<int>(round(std::max(std::max(proj_abc[0].y, proj_abc[1].y), proj_abc[2].y) - 0.5));

    double x_G = (proj_abc[0].x + proj_abc[1].x + proj_abc[2].x)/3;
    double y_G = (proj_abc[0].y + proj_abc[1].y + proj_abc[2].y)/3;
    double z_breuk_G = (1/(3*A.z)) + (1/(3*B.z)) + (1/(3*C.z));

    for(int i = min; i <= max; ++i) {
        double xl_AB = std::numeric_limits<double>::infinity();
        double xl_BC = std::numeric_limits<double>::infinity();
        double xl_AC = std::numeric_limits<double>::infinity();
        double xr_AB = - std::numeric_limits<double>::infinity();
        double xr_BC = - std::numeric_limits<double>::infinity();
        double xr_AC = - std::numeric_limits<double>::infinity();

        // Eerst voor lijn AB
        if((i - proj_abc[0].y)*(i - proj_abc[1].y) <= 0 and proj_abc[0].y != proj_abc[1].y) {
            xl_AB = xr_AB = proj_abc[1].x + ((proj_abc[0].x - proj_abc[1].x) * ((i - proj_abc[1].y)/ (proj_abc[0].y - proj_abc[1].y)));
        }
        // Eerst voor lijn AC
        if((i - proj_abc[0].y)*(i - proj_abc[2].y) <= 0 and proj_abc[0].y != proj_abc[2].y) {
            xl_AC = xr_AC = proj_abc[2].x + ((proj_abc[0].x - proj_abc[2].x) * ((i - proj_abc[2].y)/ (proj_abc[0].y - proj_abc[2].y)));
        }
        // Eerst voor lijn BC
        if((i - proj_abc[1].y)*(i - proj_abc[2].y) <= 0 and proj_abc[1].y != proj_abc[2].y) {
            xl_BC = xr_BC = proj_abc[2].x + ((proj_abc[1].x - proj_abc[2].x) * ((i - proj_abc[2].y)/ (proj_abc[1].y - proj_abc[2].y)));
        }

        int de_xl = static_cast<int> (round(std::min(std::min(xl_AB, xl_BC), xl_AC) + 0.5));
        int de_xr = static_cast<int> (round(std::max(std::max(xr_AB, xr_BC), xr_AC) - 0.5));

        for(int j = de_xl; j <= de_xr; ++j) {
            Vector3D u = B-A;
            Vector3D v = C-A;
            double w_1 = u.y * v.z - u.z * v.y;
            double w_2 = u.z * v.x - u.x * v.z;
            double w_3 = u.x * v.y - u.y * v.x;
            double k = w_1 * A.x + w_2 * A.y + w_3 * A.z;
            double dzdx = w_1/(-d * k);
            double dzdy = w_2/(-d * k);
            double de_z_breuk = 1.0001 * z_breuk_G + (j - x_G)*dzdx + (i - y_G)*dzdy;

            double z_point = 1/de_z_breuk;
            double y_point = (-z_point)*((i-dy)/d);
            double x_point = (-z_point)*((j-dx)/d);
            Vector3D total_point = Vector3D::point(x_point, y_point, z_point);
            Color temp(0,0,0);
            for (Light* light: lights) {
                // Pointlight for diffuselight
                if(auto* pointlight = dynamic_cast<PointLight*>(light)) {
                    Vector3D l = Vector3D::normalise(pointlight->location - total_point);
                    double cosalfa = (n.x * l.x) + (n.y * l.y) + (n.z * l.z);
                    double cosbeta =  cos(pointlight->spotAngle * M_PI/180);
                    if(cosalfa > cosbeta) {
                        double costheta = 1-((1-cosalfa)/ (1 - cosbeta));
                        Color diffuse = Color(diffuseReflection.red * light->diffuseLight.red *costheta, diffuseReflection.green * light->diffuseLight.green *costheta, diffuseReflection.blue * light->diffuseLight.blue *costheta);
                        temp = tooBright(temp, diffuse);
                    }
                    //Pointlight for Specular light
                    Vector3D refl_vec = Vector3D::normalise((2 * n * cosalfa) - l);
                    Vector3D pic = Vector3D::normalise(pointlight->eye_trans - total_point);
                    double costheta = Vector3D::dot(refl_vec, pic);
                    if(costheta > 0) {
                        double power = pow(costheta, reflectionCoeff);
                        Color specular = Color(power * specularReflection.red * light->specularLight.red, power * specularReflection.green * light->specularLight.green, power * specularReflection.blue * light->specularLight.blue);
                        temp = tooBright(temp, specular);
                    }
                }
                // infinity light Specular light
                if(auto* inflight = dynamic_cast<InfLight*>(light)) {
                    Vector3D l = - Vector3D::normalise(inflight->ldVector);
                    double cosalfa = (n.x * l.x) + (n.y * l.y) + (n.z * l.z);
                    Vector3D refl_vec = Vector3D::normalise((2 * n * cosalfa) - l);
                    Vector3D pic = Vector3D::normalise(inflight->eye_trans - total_point);
                    double cosbeta = Vector3D::dot(refl_vec, pic);
                    if(cosbeta > 0) {
                        double power = pow(cosbeta, reflectionCoeff);
                        Color specular = Color(power * specularReflection.red * light->specularLight.red, power * specularReflection.green * light->specularLight.green, power * specularReflection.blue * light->specularLight.blue);
                        temp = tooBright(temp, specular);
                    }
                }
            }
            if (check_smaller(de_z_breuk, i, j)) {
                temp = tooBright(temp, total_light);
                (image)(j, i) = convertcolor(temp);
            }
        }
    }
}