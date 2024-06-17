// if the Polygon class name conflicts with a class in wingdi.h on Windows, use
// a namespace or change the name
#include "Vector.h"
#include <iostream>
typedef std::pair<Vector, Vector> Edge;

class Polygon {
  public:
    std::vector<Vector> vertices;

    Polygon() {}

    Polygon(const std::vector<Vector> &vertices) : vertices(vertices) {}
    Polygon(std::vector<Vector> &&vertices) : vertices(vertices) {}
    // Polygon(std::default_random_engine &gen) { random_polygon(gen); }a

    Polygon(const Polygon &other);
    Polygon(Polygon &&other);

    // Polygon &&clip_polygon(const Polygon &clip_polygon);

    void add(const Vector &__value);

    void add(Vector &&__value);

    // std::pair<Vector, Vector> intersect(size_t i, size_t j);

    // move constructor
    // Polygon(Polygon &&other);

    void operator=(Polygon &&other);
    void operator=(const Polygon &other);

    const Vector &operator[](size_t a) const;
    Vector &operator[](size_t a);

    size_t size() const;
};

class OptimalTransportProperties {
  public:
    std::vector<Vector> &points;
    double desired_fluid_volume = 1;

    OptimalTransportProperties(std::vector<Vector> &points,
                               double desired_fluid_volume)
        : points(points), desired_fluid_volume(desired_fluid_volume) {}
};

std::vector<Vector> random_polygon(std::default_random_engine &gen, size_t n);

std::vector<Polygon> voronoy_diagram(const std::vector<Vector> &points);

bool inside_voronoy(const Vector &p, const Vector &p0, const Vector &p1);

Vector intersect_voronoy(const Vector &seg1, const Vector &seg2,
                         const Vector &p0, const Vector &p1);

Vector intersect_power(const Vector &A, const Vector &B, const Vector &p0,
                       double w0, const Vector &p1, double w1);

double *random_weights(std::default_random_engine &gen, size_t);

void voronoy_clip(Polygon &poly, const Vector &a, const Vector &b);

void power_clip(Polygon &poly, const Vector &p0, double w0, const Vector &p1,
                double w1);

std::vector<Polygon> power_diagram(const std::vector<Vector> &points,
                                   const double *weights);

double compute_area(const Polygon &polygon);
double compute_area_signed(const Polygon &polygon);

double compute_weird_formula(const Polygon &polygon, const Vector &point);

double *optimal_transport(std::vector<Vector> &points,
                          std::default_random_engine &gen);
double *optimal_transport_fluids(std::vector<Vector> &points,
                                 std::default_random_engine &gen);

void run_time_step(std::vector<Vector> &points, std::vector<Vector> &velocity,
                   double *weights, std::default_random_engine &gen,
                   const double mass, const double eps, const double dt,
                   int frame_iter);

Vector compute_centroid(const Polygon &polygon);

Polygon make_circle(const Vector &cirlce, double radius);
std::vector<Polygon> power_diagram_circle(const std::vector<Vector> &points,
                                          const double *weights);

void draw_random_polygon(std::default_random_engine &gen);
// saves a static svg file. The polygon vertices are supposed to be in the range
// [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename,
              std::string fillcol = "none");

// Adds one frame of an animated svg file. frameid is the frame number (between
// 0 and nbframes-1). polygons is a list of polygons, describing the current
// frame. The polygon vertices are supposed to be in the range [0..1], and a
// canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons,
                       std::string filename, int frameid, int nbframes);

void save_frame(const std::vector<Polygon> &cells, std::string filename,
                int frameid = 0);
int sgn(double val);

