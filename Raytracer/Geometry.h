#ifndef INTERSECT_FILE
#define INTERSECT_FILE 1
#include "Intersect.h"
#endif

#ifndef RAYTRACER_FILE
#define RAYTRACER_FILE
#include "raytracer.h"
#endif

class Geometry {
public:
    Geometry() {};
    virtual Intersect intersect(const Ray &ray) const =0;
};
