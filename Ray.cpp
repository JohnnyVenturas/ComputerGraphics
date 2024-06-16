#include "Ray.h"
#include "Scene.h"

Ray::Ray(const Camera &camera, int x, int y, bool antialiasing,
         std::default_random_engine &gen,
         std::uniform_real_distribution<double> &dist) {
    build_ray(camera, x, y, antialiasing, gen, dist);
}

Ray::Ray(const Vector &Q, const Vector &P) { build_ray(Q, P); }

void Ray::build_ray(const Camera &camera, int x, int y, bool antialiasing,
                    std::default_random_engine &gen,
                    std::uniform_real_distribution<double> &dist) {
    const Vector &Q = camera.Q;

    Vector P = Vector(Q.data[0] + x + 0.5 - camera.w / 2.,
                      Q.data[1] + y + 0.5 - camera.h / 2.,
                      Q.data[2] - camera.w / (2 * std::tan(camera.alpha / 2)));

    // antialiasing
    if (antialiasing) {
        double r1 = dist(gen), r2 = dist(gen);
        double a = cos(2 * M_PI * r1) * sqrt(1 - r2) / 2;
        double b = sin(2 * M_PI * r1) * sqrt(1 - r2) / 2;

        P += Vector(a, b, 0);
    }

    u = P - camera.Q;
    O = camera.Q;

    u.normalize();
}

void Ray::build_ray(const Vector &Q, const Vector &P) {
    // Build Ray starting point
    O = Q;
    // Give ray direction
    u = P - Q;
    u.normalize();
}

Ray Ray::build_reflected_ray(Intersect &intersection_point) const {

    // std::cout << intersection_point.N.norm() << "\n";
    Ray reflected_ray;
    Vector reflected_direction =
        u - 2 * dot(u, intersection_point.N) * intersection_point.N;

    reflected_ray.O = intersection_point.P + epsilon * intersection_point.N;
    reflected_direction.normalize();

    reflected_ray.u = reflected_direction;

    // std::cout << dot(intersection_point.N, reflected_ray.u) <<"\n";

    return reflected_ray;
}

Ray Ray::build_refracted_ray(Intersect &intersection_point) const {
    double n1 = intersection_point.object->object_properties.n1;
    double n2 = intersection_point.object->object_properties.n2;
    intersection_point.N =
        (intersection_point.object->object_properties.invert ? -1 : 1) *
        intersection_point.N;

    double __dot_product = dot(u, intersection_point.N);

    int sign = __dot_product < 0 ? -1 : 1;

    if (__dot_product > 0) {
        std::swap(n1, n2);
    }

    Vector refracted_tangent =
        (n1 / n2) * (u - __dot_product * intersection_point.N);

    double refracted_internal_angle =
        1 - (n1 / n2) * (n1 / n2) * (1 - __dot_product * __dot_product);

    Ray refracted_ray;

    if (refracted_internal_angle < 0.) {
        return build_reflected_ray(intersection_point);
    }

    Vector refracted_normal =
        (sign)*intersection_point.N * sqrt(refracted_internal_angle);

    refracted_ray.u = refracted_normal + refracted_tangent;

    refracted_ray.O =
        intersection_point.P + sign * intersection_point.N * epsilon;

    return refracted_ray;
}

void Ray::randomize_direction() {}

