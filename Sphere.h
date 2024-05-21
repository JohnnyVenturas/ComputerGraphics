#ifndef OBJECT_PROPERTIES_FILE
#define OBJECT_PROPERTIES_FILE 1
#include "ObjectProperties.h"
#endif


#ifndef INTERSECT_FILE
#define INTERSECT_FILE 1
#include "Intersect.h"
#endif

#ifndef GEOMETRY_FILE
#define GEOMTERY_FILE
#include "Geometry.h"
#endif

class Sphere : virtual public Geometry {
public:
  Vector C;
  double R;
  ObjectProperties object_properties;

  bool invert = false;

  MaterialType surface_type;

  Sphere() {}

  Sphere(Vector &C, double R) : C(C), R(R) {}

  Sphere(const Vector &C, double R, const ObjectProperties &object_properties)
      : C(C), R(R), object_properties(object_properties) {}
  Sphere(const Vector &C, double R, const ObjectProperties &object_properties,
         bool invert)
      : C(C), R(R), object_properties(object_properties), invert(invert) {}


  virtual Intersect intersect(const Ray &ray) const;
};
