
#ifndef VECTOR_FILE 
#define VECTOR_FILE 1
#include "Vector.h"
#endif

#ifndef OBJECT_PROPERTIES_FILE
#define OBJECT_PROPERTIES_FILE 1
#include "ObjectProperties.h"
#endif

#ifndef SPHERE_FILE
#define SPHERE_FILE 1
#include "Sphere.h"
#endif

class Intersect {
public:
  bool intersect;
  Vector P, N;
  ObjectProperties object_properties;
  double t;
  int objectID;
  Sphere sphere;

  Intersect() {}

  Intersect(bool intersect) : intersect(intersect) { t = 100000000; }

  Intersect(bool intersect, const Vector &P, const Vector &N,
            const ObjectProperties &object_properties, double t)
      : intersect(intersect), P(P), N(N), object_properties(object_properties),
        t(t) {}

  Intersect(bool intersect, const Vector &P, const Vector &N,
            const ObjectProperties &object_properties, double t, int id)
      : intersect(intersect), P(P), N(N), object_properties(object_properties),
        t(t), objectID(id) {}

  Intersect(bool intersect, const Vector &P, const Vector &N,
            const ObjectProperties &object_properties, double t,
            const Sphere &sphere)
      : intersect(intersect), P(P), N(N), object_properties(object_properties),
        t(t), sphere(sphere) {}
};
