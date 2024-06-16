
#ifndef SCENE_H
#define SCENE_H

#include "ObjectProperties.h"
#include "Ray.h"
#include "Vector.h"
#include <iostream>
#include <vector>

#define REFLECTION 1
#define REFRACTION 2
#define BIFRACTION 3

static std::default_random_engine gen(10);
static std::uniform_real_distribution<double> dist(0.0, 1.0);
constexpr double epsilon = 1e-7;

class Geometry {
  public:
    ObjectProperties object_properties;

    Geometry(){};
    Geometry(const ObjectProperties &object_properties)
        : object_properties(object_properties) {}

    virtual Intersect intersect(const Ray &ray) const = 0;
};

class Sphere : virtual public Geometry {
  public:
    Vector C;
    double R;

    // ObjectProperties object_properties;
    // MaterialType surface_type;

    Sphere() {}

    Sphere(Vector &C, double R) : C(C), R(R) {}

    Sphere(const Vector &C, double R, const ObjectProperties &object_properties)
        : Geometry(object_properties), C(C), R(R) {}

    virtual Intersect intersect(const Ray &ray) const;
};

class LightSource {
  public:
    // Light Source
    Vector center;
    // Intensity;
    double intensity;

    LightSource() {}

    LightSource(const Vector &light_source, double intensity)
        : center(light_source), intensity(intensity) {}
};

class Scene {
  public:
    std::vector<Geometry *> scene;
    LightSource light_source;

    Scene() {}

    Scene(const std::vector<Geometry *> &objects) : scene(objects) {}

    Scene(const std::vector<Geometry *> &objects,
          const LightSource light_source)
        : scene(objects), light_source(light_source) {}

    // Get color from interesection with light and object

    // outside is true if the interesecting ray is emitted from the outside of
    // the object
    Vector get_color(const Ray &ray, int depth);

    Intersect intersect(const Ray &ray) const;
};



#endif
