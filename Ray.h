#include "Camera.h"
#include "Intersect.h"
#include "Vector.h"
#ifndef RAY_H
#define RAY_H
class Ray {
  public:
    Vector O, u;
    Ray() {}

    // Standard way to build a ray, from an origin_point and a direction
    // Ray(const Vector &O, const Vector &u) : O(O), u(u) {}

    // Build a ray from a camera, and a coordinate on the camera screen
    Ray(const Camera &camera, int x, int y, bool antialiasing,
        std::default_random_engine &engine,
        std::uniform_real_distribution<double> &dist);

    // Build a ray from an origin point and end point
    Ray(const Vector &Q, const Vector &P);

    // Build a ray from an intersection_point and some refractive/reflexive/etc
    // surface
    Ray(const Ray &ray, const Intersect &intersection_point);

    // Build a ray from a camera, and a coordinates on the camera screen
    void build_ray(const Camera &camera, int x, int y, bool antialiasing,
                   std::default_random_engine &engine,
                   std::uniform_real_distribution<double> &dist);

    // Build a ray from an origin point and end point
    void build_ray(const Vector &Q, const Vector &P);

    Ray build_refracted_ray(Intersect &intersection_point) const;
    Ray build_reflected_ray(Intersect &intersection_point) const;

    void randomize_direction();

    // Build a ray from an intersection_point and some refractive/reflexive/etc
    // surface

    // Polymorphic interesection function of a ray
    // Implemented for spheres and scene
    // template <typename T> Intersect get_intersection(const T &object) const;
};

#endif
