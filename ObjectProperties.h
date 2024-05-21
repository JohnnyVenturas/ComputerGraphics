
#ifndef VECTOR_FILE
#define VECTOR_FILE 1
#include "Vector.h"
#endif
enum MaterialType {
  DIFUSIVE = 0,
  REFLECTIVE = 1,
  REFRACTIVE = 2,
  BIFRACTIVE = 3,
};

class ObjectProperties {
public:
  Vector albedo;
  MaterialType material_type;

  double n1 = 1, n2 = 1.5; // n1 refractive surface for outside the object, n2
                           // for inside the ojbect

  ObjectProperties() {}
  explicit ObjectProperties(const Vector &albedo)
      : albedo(albedo), material_type(DIFUSIVE) {}
  ObjectProperties(const Vector &albedo, const MaterialType material_type)
      : albedo(albedo), material_type(material_type) {}

  ObjectProperties(const Vector &albedo, const MaterialType material_type,
                   double n1, double n2)
      : albedo(albedo), material_type(material_type), n1(n1), n2(n2) {}
};
