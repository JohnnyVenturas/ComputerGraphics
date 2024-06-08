#include <cmath>
#include <random>
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION

#ifndef RAYTRACER_FILE
#define RAYTRACER_FILE 1
#include "raytracer.h"
#endif

#include "stb_image.h"
#include <iostream>
#define RAY_NUMBER 128

#define RAYTRACING 1
#define FRESNEL 1
#define ANTIALIASING 1
#define CAT 0
#define INTERPOLATION 1
#define LIGHT_BOUNCE 6

#include "mesh.h"

// static std::random_device rd;
static std::default_random_engine gen(10);
static std::uniform_real_distribution<double> dist(0.0, 1.0);
double epsilon = 1e-7;
Vector::Vector(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

Vector Vector ::sample_direction() {

    Vector normal = *this;
    double r1 = dist(gen), r2 = dist(gen);

    Vector T1 = normal;

    int pos = 0;

    if (abs(T1.data[0]) <= abs(T1.data[1]) &&
        abs(T1.data[0]) <= abs(T1.data[2])) {
        pos = 0;
    } else if (abs(T1.data[1]) <= abs(T1.data[2]) &&
               abs(T1.data[1]) <= abs(T1.data[0])) {
        pos = 1;
    } else {
        pos = 2;
    }
    T1.data[pos] = 0;

    std::swap(T1.data[(pos + 1) % 3], T1.data[(pos + 2) % 3]);

    T1.normalize();

    T1.data[(pos + 1) % 3] = -T1.data[(pos + 1) % 3];

    Vector T2 = cross(normal, T1);

    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);
    Vector res = x * T1 + y * T2 + z * normal;
    return res;
}

void Vector::operator+=(const Vector &other) {
    this->data[0] += other.data[0];
    this->data[1] += other.data[1];
    this->data[2] += other.data[2];
}
void print(const Vector &a, const Vector &b) {
    printf("%f %f %f\n", a.data[0], a.data[1], a.data[2]);
    printf("%f %f %f\n", b.data[0], b.data[1], b.data[2]);
}

void print(const Vector &a) {
    printf("%f %f %f\n", a.data[0], a.data[1], a.data[2]);
}

double Vector::norm2() const {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
}

double Vector::norm() const { return sqrt(norm2()); }
void Vector::normalize() {
    double n = norm();
    data[0] /= n;
    data[1] /= n;
    data[2] /= n;
}

double Vector ::operator[](int i) const { return data[i]; }
double &Vector::operator[](int i) { return data[i]; }
bool Vector::operator==(const Vector &other) const {
    return abs(data[0] - other.data[0]) < epsilon &&
           abs(data[1] - other.data[1]) < epsilon &&
           abs(data[2] - other.data[2]) < epsilon;
}

bool Vector::operator!=(const Vector &other) const { return !(*this == other); }

Vector Vector::operator-() const {
    return Vector(-data[0], -data[1], -data[2]);
}

Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator*(const Vector &a, const Vector &b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const double a, const Vector &b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator/(const Vector &a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

double dot(const Vector &a, const Vector &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector &a, const Vector &b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
}

bool operator>=(const Vector &a, const Vector &b) {
    if (a == b)
        return true;
    if (a[0] < b[0])
        return false;
    if (a[1] < b[1])
        return false;
    if (a[2] < b[2])
        return false;
    return true;
}

bool operator<=(const Vector &a, const Vector &b) {
    if (a == b)
        return true;
    return !(a >= b);
}

Ray::Ray(const Camera &camera, int x, int y, bool antialiasing) {
    build_ray(camera, x, y, antialiasing);
}
Ray::Ray(const Vector &Q, const Vector &P) { build_ray(Q, P); }

void Ray::build_ray(const Camera &camera, int x, int y, bool antialiasing) {
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

     //std::cout << "therre are \n";
     //std::cout << intersect_point.object << " ";
     //std::cout <<"muie la cur\n";

    return intersect_point;
}

Ray Ray::build_reflected_ray(Intersect &intersection_point) const {

    //std::cout << intersection_point.N.norm() << "\n";
    Ray reflected_ray;
    Vector reflected_direction =
        u - 2 * dot(u, intersection_point.N) * intersection_point.N;

    reflected_ray.O = intersection_point.P + epsilon * intersection_point.N;
    reflected_direction.normalize();

    reflected_ray.u = reflected_direction;

    //std::cout << dot(intersection_point.N, reflected_ray.u) <<"\n";


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

Vector Scene::get_color(const Ray &ray, int depth) {


    if (depth < 0)
        return Vector(0, 0, 0);

    Intersect intersection_point = this->intersect(ray);
    if(!intersection_point.intersect) return Vector(0,0,0);

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
        double check = scene_intersect.t + epsilon< distance_to_light ? 0. : 1.;

        Vector L0 = light_source.intensity / (4 * M_PI * d * d) * 1 / M_PI *
               check *
                      std::max(dot(intersection_point.N, -light_ray.u), 0.) *
                      intersection_point.object->object_properties.albedo;

        Vector ray_end = intersection_point.N.sample_direction() +
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

void draw_image(int W, int H) {
    Camera camera = Camera(Vector(0, 0, 55), W, H, M_PI / 3);
    LightSource light_source(Vector(0, 0, 20), 1e10 / 4);
    Vector color_white = Vector(1, 1, 1);
    Vector color_grey = Vector(.5, .5, .5);
    Vector color_green(0, 128. / 255, 0), color_magenta(1, 0, 1),
        color_blue(0, 0, 1), color_red(1, 0, 0), color_yellow(1, 1, 0),
        color_cyan(0, 1, 1), color_cat(0.3, 0.2, 0.25);

    std::vector<Geometry *> objects;

#if FRESNEL == 1
    ObjectProperties reflective = ObjectProperties(color_white, REFLECTIVE);
    ObjectProperties hollow = ObjectProperties(color_white, BIFRACTIVE, 1, 1.5);
    ObjectProperties bifractive = ObjectProperties(color_white, BIFRACTIVE);
#else
    ObjectProperties reflective = ObjectProperties(color_white, REFLECTIVE);
    ObjectProperties hollow = ObjectProperties(color_white, REFRACTIVE, 1, 1.5);
    ObjectProperties bifractive = ObjectProperties(color_white, REFRACTIVE);
#endif

    ObjectProperties difusive_green = ObjectProperties(color_green, DIFUSIVE);
    ObjectProperties difusive_magenta =
        ObjectProperties(color_magenta, DIFUSIVE);
    ObjectProperties difusive_blue = ObjectProperties(color_blue, DIFUSIVE);
    ObjectProperties difusive_red = ObjectProperties(color_red, DIFUSIVE);
    ObjectProperties difusive_yellow = ObjectProperties(color_yellow, DIFUSIVE);
    ObjectProperties difusive_cyan = ObjectProperties(color_cyan, DIFUSIVE);
    ObjectProperties difusive_white = ObjectProperties(color_white, DIFUSIVE);
    ObjectProperties difusive_grey = ObjectProperties(color_grey, DIFUSIVE);
    ObjectProperties difusive_cat = ObjectProperties(color_cat, DIFUSIVE);

    Sphere reflective_sphere(Vector(-20, 0, 0), 10, reflective);
    Sphere hollow_sphere(Vector(20, 0, 0), 10, hollow);

    // inner object(should not have predefined prop)
    Sphere hollow_inner(
        Vector(20, 0, 0), 9.5,
        ObjectProperties(color_white, REFRACTIVE, 1, 1.5, true));

    // target sphere

    Sphere refractive_sphere(Vector(0, 0, 0), 10, bifractive);

    Sphere wall_front(Vector(0, 0, -1000), 940, difusive_green);

    Sphere wall_behind(Vector(0, 0, 1000), 940, difusive_magenta);

    Sphere wall_bottom(Vector(0, -1000, 0), 990, difusive_blue);

    Sphere wall_top(Vector(0, 1000, 0), 940, difusive_red);

    Sphere wall_right(Vector(1000, 0, 0), 940, difusive_yellow);

    Sphere wall_left(Vector(-1000, 0, 0), 940, difusive_cyan);
#if CAT == 0
    objects.push_back(&reflective_sphere);
    objects.push_back(&hollow_sphere);
    objects.push_back(&hollow_inner);
    objects.push_back(&refractive_sphere);
#endif
    objects.push_back(&wall_left);
    objects.push_back(&wall_right);
    objects.push_back(&wall_top);
    objects.push_back(&wall_bottom);
    objects.push_back(&wall_front);
    objects.push_back(&wall_behind);

#if CAT == 1
    TriangleMesh *cat_mesh = new TriangleMesh("./cat/cat.obj", reflective);

    cat_mesh ->transform(0.6, Vector(0, -10, -10));
    cat_mesh -> root = new BVH(*cat_mesh);
    objects.push_back(cat_mesh);
#endif

    Scene scene(objects, light_source);

    std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector color(0, 0, 0);

            for (int k = 0; k < RAY_NUMBER; ++k) {
                int x = j, y = H - i - 1;

                // true enables anti-aliasing, revert to false otherwise
                Ray camera_ray(camera, x, y, ANTIALIASING);

                color += scene.get_color(camera_ray, LIGHT_BOUNCE - 1);
            }

            color = color / RAY_NUMBER;
            double gamma = 1 / 2.2;

            image[(i * W + j) * 3 + 0] =
                std::min(std::pow(color.data[0], gamma), 255.0);
            image[(i * W + j) * 3 + 1] =
                std::min(std::pow(color.data[1], gamma), 255.0);
            image[(i * W + j) * 3 + 2] =
                std::min(std::pow(color.data[2], gamma), 255.0);
        }   
    }

     //std::cout <<"Sanity check\n";
     //std::cout << static_cast<Geometry *>(cat_mesh) << " " << cat_mesh -> object_properties.material_type<< "\n";
     //std::cout <<"End\n";

    stbi_write_png("image.png", W, H, 3, &image[0], 0);
}

int main() {
    int W = 512;
    int H = 512;

    draw_image(W, H);

    return 0;
}
