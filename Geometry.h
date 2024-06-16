#ifndef GEOMETRY_H
#define GEOMETRY_H 
#include "ObjectProperties.h"
#include "Intersect.h"
#include "Ray.h"
class Geometry {
  public:
    ObjectProperties object_properties;

    Geometry(){};
    Geometry(const ObjectProperties &object_properties)
        : object_properties(object_properties) {}

    virtual Intersect intersect(const Ray &ray) const = 0;
};

#endif
