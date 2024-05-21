#include <iostream>
#include <vector>

#define REFLECTION 1
#define REFRACTION 2
#define BIFRACTION 3

class Vector {
  public:
    explicit Vector(double x = 0, double y = 0, double z = 0);
    double norm2() const;
    double norm() const;
    void normalize();

    double operator[](int i) const;
    double &operator[](int i);

    Vector operator-() const;
    void operator+=(const Vector &other);

    bool operator==(const Vector &other) const;
    bool operator!=(const Vector &other) const;
    Vector sample_direction();

    // void operator=(const Vector &other) const;

    double data[3];
};

Vector cross(const Vector &a, const Vector &b);
Vector operator+(const Vector &a, const Vector &b);
Vector operator*(const Vector &a, const Vector &b);
Vector operator-(const Vector &a, const Vector &b);
Vector operator*(const double a, const Vector &b);
Vector operator*(const Vector &a, const double b);
Vector operator/(const Vector &a, const double b);

bool operator>=(const Vector &a, const Vector &b);
bool operator<=(const Vector &a, const Vector &b);

double dot(const Vector &a, const Vector &b);
Vector cross(const Vector &a, const Vector &b);

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
    bool invert = false;

    ObjectProperties() {}
    explicit ObjectProperties(const Vector &albedo)
        : albedo(albedo), material_type(DIFUSIVE) {}
    ObjectProperties(const Vector &albedo, const MaterialType material_type)
        : albedo(albedo), material_type(material_type) {}

    ObjectProperties(const Vector &albedo, const MaterialType material_type,
                     double n1, double n2)
        : albedo(albedo), material_type(material_type), n1(n1), n2(n2) {}

    ObjectProperties(const Vector &albedo, const MaterialType material_type,
                     double n1, double n2, bool invert)
        : albedo(albedo), material_type(material_type), n1(n1), n2(n2),
          invert(invert) {}
};

class Intersect;
class Ray;

class Geometry {
  public:
    ObjectProperties object_properties;

    Geometry(){};
    Geometry(const ObjectProperties &object_properties)
        : object_properties(object_properties) {}

    virtual Intersect intersect(const Ray &ray) const = 0;
};

class Camera {
  public:
    Vector Q;
    int w, h;
    double alpha;

    Camera(const Vector &Q, int w, int h, double alpha)
        : Q(Q), w(w), h(h), alpha(alpha) {}
};

class Ray {
  public:
    Vector O, u;
    Ray() {}

    // Standard way to build a ray, from an origin_point and a direction
    // Ray(const Vector &O, const Vector &u) : O(O), u(u) {}

    // Build a ray from a camera, and a coordinate on the camera screen
    Ray(const Camera &camera, int x, int y, bool antialiasing);

    // Build a ray from an origin point and end point
    Ray(const Vector &Q, const Vector &P);

    // Build a ray from an intersection_point and some refractive/reflexive/etc
    // surface
    Ray(const Ray &ray, const Intersect &intersection_point);

    // Build a ray from a camera, and a coordinates on the camera screen
    void build_ray(const Camera &camera, int x, int y, bool antialiasing);

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

class Sphere : virtual public Geometry {
  public:
    Vector C;
    double R;

    // ObjectProperties object_properties;
    // MaterialType surface_type;

    Sphere() {}

    Sphere(Vector &C, double R) : C(C), R(R) {}

    Sphere(const Vector &C, double R, const ObjectProperties &object_properties)
        : Geometry(object_properties), C(C), R(R) {}

    virtual Intersect intersect(const Ray &ray) const;
};

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

class LightSource {
  public:
    // Light Source
    Vector center;
    // Intensity;
    double intensity;

    LightSource() {}

    LightSource(const Vector &light_source, double intensity)
        : center(light_source), intensity(intensity) {}
};

class Scene {
  public:
    std::vector<Geometry *> scene;
    LightSource light_source;

    Scene() {}

    Scene(const std::vector<Geometry *> &objects) : scene(objects) {}

    Scene(const std::vector<Geometry *> &objects,
          const LightSource light_source)
        : scene(objects), light_source(light_source) {}

    // Get color from interesection with light and object

    // outside is true if the interesecting ray is emitted from the outside of
    // the object
    Vector get_color(const Ray &ray, int depth);

    Intersect intersect(const Ray &ray) const;
};
