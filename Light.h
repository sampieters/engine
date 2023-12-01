//
// Created by Sam Pieters on 06/05/2020.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H

#include <list>
#include "Color.h"
#include "vector3d.h"
#include "easy_image.h"

class Light {
public:
    //de ambiente licht component
    Color ambientLight;
    //de diffuse licht component
    Color diffuseLight;
    //de specular licht component
    Color specularLight;
    // De eye trans voor specular
    Vector3D eye_trans;

    Light(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight);

    virtual ~Light();
};

typedef std::list<Light*> Lights3D;

class ZBuffer: public std::vector<std::vector <double>> {
public:
    std::vector<std::vector<double>> buffer;
//Constructor: maakt een Z-Buffer van de correcte
//grootte aan en initialiseert alle velden op +inf
    ZBuffer(const int width, const int height);

    void draw_zbuf_line(img::EasyImage & image, const unsigned int x0, const unsigned int y0, const double z0,
                        const unsigned int x1, const unsigned int y1, const double z1, const Color &color);

    double z_vergelijking(double z_a, double z_b, int i, int a);

    bool check_smaller(double z, int row, int column);

    void draw_zbuf_triag(img::EasyImage& image, Vector3D const& A, Vector3D const& B, Vector3D const& C, double d, double dx, double dy, Color ambientReflection, Color diffuseReflection, Color specularReflection, double reflectionCoeff, Lights3D& lights);
};

class InfLight: public Light {
public:
    //de richting waarin het licht schijnt
    Vector3D ldVector;

    InfLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
             const Vector3D &ldVector);
};

class PointLight: public Light {
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;

    //ZBuffer shadowMask;

    //De eyepoint transformatie matrix om de wereld vanuit het oogpunt van het lichtpunt weer te geven 
    //Matrix eye;
    //Waardes nodig om een punt in de shadowmask op te zoeken
    //double d, dx, dy;

    PointLight(const Color &ambientLight, const Color &diffuseLight, const Color &specularLight,
               const Vector3D &location, double spotAngle);
};





#endif //ENGINE_LIGHT_H
