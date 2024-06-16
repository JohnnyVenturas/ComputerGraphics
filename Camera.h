#include "Vector.h"
#ifndef CAMERA_H
#define CAMERA_H
class Camera {
  public:
    Vector Q;
    int w, h;
    double alpha;

    Camera(const Vector &Q, int w, int h, double alpha)
        : Q(Q), w(w), h(h), alpha(alpha) {}
};

#endif
