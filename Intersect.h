#ifndef INTERSECT_H
#define INTERSECT_H
#include "Vector.h"

class Geometry;
class Intersect {
  public:
    bool intersect;
    Vector P, N;
    //ObjectProperties object_properties;
    double t;
    int objectID;
    const Geometry *object;

    Intersect() {}

    Intersect(bool intersect) : intersect(intersect) { t = 100000000; }

    Intersect(bool intersect, const Vector &P, const Vector &N,
              double t)
        : intersect(intersect), P(P), N(N),
          t(t) {}

    Intersect(bool intersect, const Vector &P, const Vector &N,
               double t, int id)
        : intersect(intersect), P(P), N(N),
           t(t), objectID(id) {}

    Intersect(bool intersect, const Vector &P, const Vector &N,
              double t,
              const Geometry *object)
        : intersect(intersect), P(P), N(N),
           t(t), object(object) {}
};

#endif
