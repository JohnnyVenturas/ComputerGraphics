#include "Scene.h"
#include "Intersect.h"
#include "raytracer.h"
Intersect Sphere::intersect(const Ray &ray) const {
    double __left_term = dot(ray.u, C - ray.O);
    double __delta_dot = dot(ray.u, ray.O - C);
    double delta = __delta_dot * __delta_dot - ((ray.O - C).norm2() - R * R);
    if (delta < 0)
        return Intersect(false);

    double t1 = __left_term - sqrt(delta);
    double t2 = __left_term + sqrt(delta);

    if (t1 < 0 && t2 < 0)
        return Intersect(false);

    if (t1 < 0)
        std::swap(t1, t2);

    double t = t1;
    Vector P = ray.O + t * ray.u;
    Vector N = (P - C);
    N.normalize();
    return Intersect(true, P, N, t, static_cast<const Geometry *>(this));
}

Intersect Scene::intersect(const Ray &ray) const {
    Intersect intersect_point(false);
    for (auto &sphere : scene) {
        const Intersect &current_intersect = sphere->intersect(ray);
        if (!current_intersect.intersect) {
            continue;
        }

        if (current_intersect.t < intersect_point.t) {
            intersect_point = current_intersect;
        }
    }

    // std::cout << "therre are \n";
    // std::cout << intersect_point.object << " ";
    // std::cout <<"muie la cur\n";

    return intersect_point;
}

Vector Scene::get_color(const Ray &ray, int depth) {

    if (depth < 0)
        return Vector(0, 0, 0);

    Intersect intersection_point = this->intersect(ray);
    if (!intersection_point.intersect)
        return Vector(0, 0, 0);

    int material_type =
        intersection_point.object->object_properties.material_type;

    double __dot_product = dot(intersection_point.N, ray.u);

    if (material_type == BIFRACTIVE) {

        double n1 = intersection_point.object->object_properties.n1;
        double n2 = intersection_point.object->object_properties.n2;
        double k0 = ((n1 - n2) * (n1 - n2)) / ((n1 + n2) * (n1 + n2));
        double dot_product = 1 - abs(dot(intersection_point.N, ray.u));
        double dot_product2 = dot_product * dot_product;
        double dot_product4 = dot_product2 * dot_product2;
        double dot_product5 = dot_product * dot_product4;
        double R = k0 + (1 - k0) * dot_product5;
        double next = dist(gen);

        if (__dot_product > 0) {
            material_type = REFRACTIVE;
        } else {
            material_type = next < R ? REFLECTIVE : REFRACTIVE;
        }
    }

    if (material_type == DIFUSIVE) {
        double d = (intersection_point.P - light_source.center).norm();
        Ray light_ray(light_source.center,
                      intersection_point.P + intersection_point.N * epsilon);

        const Geometry *current_object = intersection_point.object;
        const Intersect scene_intersect = this->intersect(light_ray);

        double distance_to_light = (intersection_point.P - light_ray.O).norm();
        double check =
            scene_intersect.t + epsilon < distance_to_light ? 0. : 1.;

        Vector L0 = light_source.intensity / (4 * M_PI * d * d) * 1 / M_PI *
                    check *
                    std::max(dot(intersection_point.N, -light_ray.u), 0.) *
                    intersection_point.object->object_properties.albedo;

        Vector ray_end = intersection_point.N.sample_direction(gen, dist) +
                         intersection_point.P + intersection_point.N * epsilon;

        Ray random_ray =
            Ray(intersection_point.P + intersection_point.N * epsilon, ray_end);

#if RAYTRACING == 0
        return L0;
#else
        return L0 + intersection_point.object->object_properties.albedo *
                        get_color(random_ray, depth - 1);
#endif
    }

    // BUILDING REFLECTIVE RAY
    if (material_type == REFLECTIVE) {

        Ray reflected_ray = ray.build_reflected_ray(intersection_point);

        return get_color(reflected_ray, depth - 1);
    }

    // BUILDING REFRACTIVE R1000
    if (material_type == REFRACTIVE) {
        Ray refracted_ray = ray.build_refracted_ray(intersection_point);
        return get_color(refracted_ray, depth - 1);
    }

    return Vector(0, 0, 0);
}
