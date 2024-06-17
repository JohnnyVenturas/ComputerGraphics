// if the Polygon class name conflicts with a class in wingdi.h on Windows, use
// a namespace or change the name
#include "Vector.h"
#include "lbfgs.h"
#include "fluids.h"
#include <cassert>
#include <cmath>
#include <math.h>
#include <random>
#include <sstream>
#include <unistd.h>
#include <vector>
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include <iostream>

// saves a static svg file. The polygon vertices are supposed to be in the range
// [0..1], and a canvas of size 1000x1000 is created

// This might cause a bug

Polygon::Polygon(Polygon &&other) { vertices = std::move(other.vertices); }

void Polygon::operator=(Polygon &&other) {
    if (this == &other)
        return;
    vertices = std::move(other.vertices);
}

void Polygon::operator=(const Polygon &other) {
    if (this == &other)
        return;
    vertices = other.vertices;
}

Polygon::Polygon(const Polygon &other) {
    if (this == &other)
        return;
    vertices = other.vertices;
}

std::vector<Polygon> voronoy_diagram(const std::vector<Vector> &points) {
    Polygon poly(
        {Vector(0, 0, 0), Vector(0, 1, 0), Vector(1, 1, 0), Vector(1, 0, 0)});

    std::vector<Polygon> res(points.size());

#pragma omp parallel for
    for (size_t i = 0; i < points.size(); ++i) {

        Polygon out_polygon = poly;

        for (size_t j = 0; j < points.size(); ++j) {

            if (i == j) {
                continue;
            }

            voronoy_clip(out_polygon, points[i], points[j]);
        }

        res[i] = out_polygon;
    }

    return res;
}

void voronoy_clip(Polygon &poly, const Vector &p0, const Vector &p1) {

    Polygon out_polygon;

    for (size_t i = 0; i < poly.size(); ++i) {

        const Vector &cur_vertex = poly[i];
        const Vector &prev_vertex = poly[i > 0 ? (i - 1) : poly.size() - 1];

        const Vector &intersection =
            intersect_voronoy(prev_vertex, cur_vertex, p0, p1);

        if (inside_voronoy(cur_vertex, p0, p1)) {

            if (!inside_voronoy(prev_vertex, p0, p1)) {
                out_polygon.add(intersection);
            }

            out_polygon.add(cur_vertex);

        } else if (inside_voronoy(prev_vertex, p0, p1)) {

            out_polygon.add(intersection);
        }
    }

    poly = out_polygon;
}

Vector intersect_voronoy(const Vector &A, const Vector &B, const Vector &p0,
                         const Vector &p1) {
    const Vector M = (p0 + p1) / 2;

    double t = dot(M - A, p1 - p0) / dot(B - A, p1 - p0);

    return A + t * (B - A);
}

bool inside_voronoy(const Vector &p, const Vector &p0, const Vector &p1) {

    return (p - p0).norm2() <= (p - p1).norm2();
}

bool inside_power(const Vector &p, const Vector &p0, double w0,
                  const Vector &p1, double w1) {

    return ((p - p0).norm2() - w0) <= ((p - p1).norm2() - w1);
}

Vector intersect_power(const Vector &A, const Vector &B, const Vector &p0,
                       double w0, const Vector &p1, double w1) {
    // std::cout << w0 << " " << w1 << "\n";

    const Vector M = (p0 + p1) / 2;
    const Vector M1 = M + (w0 - w1) / (2 * (p0 - p1).norm2()) * (p1 - p0);

    double t = dot(M1 - A, p1 - p0) / dot(B - A, p1 - p0);

    return A + t * (B - A);
}

void power_clip(Polygon &poly, const Vector &p0, double w0, const Vector &p1,
                double w1) {

    Polygon out_polygon;
    // assert(out_polygon.size() == 0);

    for (size_t i = 0; i < poly.size(); ++i) {

        const Vector &cur_vertex = poly[i];
        const Vector &prev_vertex = poly[i > 0 ? (i - 1) : poly.size() - 1];

        const Vector &intersection =
            intersect_power(prev_vertex, cur_vertex, p0, w0, p1, w1);

        if (inside_power(cur_vertex, p0, w0, p1, w1)) {

            if (!inside_power(prev_vertex, p0, w0, p1, w1)) {
                out_polygon.add(intersection);
            }

            out_polygon.add(cur_vertex);

        } else if (inside_power(prev_vertex, p0, w0, p1, w1)) {

            out_polygon.add(intersection);
        }
    }

    poly = out_polygon;
}

std::vector<Polygon> power_diagram(const std::vector<Vector> &points,
                                   const double *weights) {

    std::vector<Polygon> res(points.size());

    Polygon poly(
        {Vector(0, 0, 0), Vector(1, 0, 0), Vector(1, 1, 0), Vector(0, 1, 0)});

#pragma omp parallel for
    for (size_t i = 0; i < points.size(); ++i) {

        Polygon out_polygon = poly;

        for (size_t j = 0; j < points.size(); ++j) {

            if (i == j) {
                continue;
            }

            power_clip(out_polygon, points[i], weights[i], points[j],
                       weights[j]);
        }

        res[i] = out_polygon;
    }

    return res;
}

std::vector<Polygon> power_diagram_circle(const std::vector<Vector> &points,
                                          const double *weights) {

    std::vector<Polygon> res(points.size());
    const size_t N = points.size();

    Polygon poly(
        {Vector(0, 0, 0), Vector(0, 1, 0), Vector(1, 1, 0), Vector(1, 0, 0)});

#pragma omp parallel for
    for (size_t i = 0; i < points.size(); ++i) {

        // std::cout << weights[N] << "\n";
        double radius = std::sqrt(weights[i] - weights[N]);

        Polygon out_polygon = make_circle(points[i], radius);

        for (size_t j = 0; j < points.size(); ++j) {

            if (i == j) {
                continue;
            }

            power_clip(out_polygon, points[i], weights[i], points[j],
                       weights[j]);
        }

        res[i] = out_polygon;
    }

    return res;
}

std::vector<Vector> random_polygon(std::default_random_engine &gen, size_t n) {

    std::uniform_real_distribution<double> dist(0, 1.);
    std::vector<Vector> vertices;

    for (int i = 0; i < n; ++i) {
        double x = 0.5 + (dist(gen) - 0.5) * 0.45;
        double y = 0.5 + (dist(gen) - 0.5) * 0.45;
        double z = 0;
        vertices.push_back(Vector(x, y, z));
    }

    return vertices;
}

double *random_weights(std::default_random_engine &gen, size_t n) {

    std::uniform_real_distribution<double> dist(0, 1.);
    /// std::vector<double> weights;

    double *weights = static_cast<double *>(malloc(n * sizeof(double)));

    for (int i = 0; i < n; ++i) {

        double x = dist(gen);
        weights[i] = 1;
    }

    weights[n - 1] = 0;

    return weights;
}

void draw_random_polygon(std::default_random_engine &gen) {
    // std::default_random_engine muie(100);

    std::vector<Vector> points = random_polygon(gen, 700);
    std::vector<Vector> velocity(points.size());

    // double *weights = optimal_transport_fluids(points, gen);
    double *weights = random_weights(gen, points.size() + 1);
    // run_time_step(points, velocity, gen, 1);

    // std::vector<Polygon> polygon_voronoy = voronoy_diagram(points);
    // save_svg(polygon_voronoy, "../tests/voronoy.svg");

    // double *weights = optimal_transport(points, gen);
    // std::vector<Polygon> polygon_power = power_diagram(points, weights);
    //
    // save_svg(polygon_power, "../tests/power.svg");

    for (int i = 0; i < 500; ++i) {
        run_time_step(points, velocity, weights, gen, 200, 2e-3, 4e-3, i);
    }

    // save_svg(polygons, "out.svg");
}

lbfgsfloatval_t evaluate(void *__points, const double *weights,
                         lbfgsfloatval_t *g, const int __number_points,
                         const lbfgsfloatval_t step) {

    const std::vector<Vector> &points =
        *static_cast<const std::vector<Vector> *>(__points);

    // assert(points.size() == __number_points);
    // assert(points.size() != 0);

    std::vector<Polygon> polygons = power_diagram(points, weights);

    double ret = 0.;
    double area = 0.;

    for (size_t i = 0; i < __number_points; ++i) {
        area = compute_area(polygons[i]);

        g[i] = area - 1. / __number_points;

        area = -area;

        ret += compute_weird_formula(polygons[i], points[i]) +
               1. / __number_points * weights[i] + weights[i] * area;
    }

    return -ret;
}

double compute_area(const Polygon &polygon) {
    double area = 0.;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const Vector &cur_vertex = polygon[i];
        const Vector &prev_vertex =
            i > 0 ? polygon[i - 1] : polygon[polygon.size() - 1];

        area +=
            (prev_vertex[0] * cur_vertex[1] - cur_vertex[0] * prev_vertex[1]);
    }
    // Test compute_area with known polygons

    return std::abs(area) / 2.;
}

double compute_weird_formula(const Polygon &polygon, const Vector &p) {
    double ret = 0.;

    for (size_t i = 0; i < polygon.size(); ++i) {
        const Vector &cur_vertex = polygon[i];
        const Vector &next_vertex =
            i < (polygon.size() - 1) ? polygon[i + 1] : polygon[0];

        double first =
            cur_vertex[0] * next_vertex[1] - next_vertex[0] * cur_vertex[1];

        double second =
            cur_vertex[0] * cur_vertex[0] + cur_vertex[0] * next_vertex[0] +
            next_vertex[0] * next_vertex[0] + cur_vertex[1] * cur_vertex[1] +
            cur_vertex[1] * next_vertex[1] + next_vertex[1] * next_vertex[1];

        double third = -4. * (p[0] * (cur_vertex[0] + next_vertex[0]) +
                              p[1] * (cur_vertex[1] + next_vertex[1]));

        double fourth = 6. * p.norm2();

        ret += first * (second + third + fourth);
    }
    return ret / 12.;
}

double *optimal_transport(std::vector<Vector> &points,
                          std::default_random_engine &gen) {

    size_t n = points.size();

    double *weights = random_weights(gen, n);

    // double fx;

    int ret = lbfgs(n, weights, NULL, evaluate, NULL,
                    static_cast<void *>(&points), NULL);

    if (ret != 0) {
        printf("MUie la cur %d\n", ret);
    }

    return weights;
}

lbfgsfloatval_t evaluate_fluids(void *instance, const double *weights,
                                lbfgsfloatval_t *g, const int __number_points,
                                const lbfgsfloatval_t step) {

    OptimalTransportProperties &properties =
        *static_cast<OptimalTransportProperties *>(instance);

    // assert(__number_points == (properties.points.size() + 1));

    const double &desired_fluid_volume = properties.desired_fluid_volume;
    std::vector<Vector> &points = properties.points;
    std::vector<Polygon> polygons = power_diagram_circle(points, weights);

    // assert(polygons.size() > 0);

    const size_t N = __number_points - 1;

    double w_air = weights[N];

    // assert(desired_fluid_volume == 0.25);

    const double desired_air_volume = 1 - desired_fluid_volume;

    double __estimated_air_volume = 1.;

    for (int i = 0; i < N; ++i) {
        __estimated_air_volume -= compute_area(polygons[i]);
    }

    // std::cout << __estimated_air_volume << "\n";
    //  assert(__estimated_air_volume >= 0);

    const double estimated_air_volume = __estimated_air_volume;

    double area = .0, ret = w_air * (desired_air_volume - estimated_air_volume);

    for (int i = 0; i < N; ++i) {
        // std::cout << weights[i] << " ";

        area = compute_area(polygons[i]);

        ret +=
            compute_weird_formula(polygons[i], points[i]) - weights[i] * area;

        ret += desired_fluid_volume / N * weights[i];

        g[i] = area - desired_fluid_volume / N;
        // std::cout << ret << "\n";
    }

    // std::cout << "\n";

    g[N] = -desired_air_volume + estimated_air_volume;

    return -ret;
}

double *optimal_transport_fluids(std::vector<Vector> &points,
                                 std::default_random_engine &gen,
                                 double *weights) {
    size_t n = points.size();

    // double *weights = random_weights(gen, n + 1);

    // for (int i = 0; i < n + 1; ++i) {
    // std::cout << weights[i] << " ";
    //}
    // std::cout << "\n";

    // double fx;

    OptimalTransportProperties properties(points, 0.25);
    int ret = lbfgs(n + 1, weights, NULL, evaluate_fluids, NULL,
                    static_cast<void *>(&properties), NULL);

    if (ret != 0) {
        printf("MUie la cur %d\n", ret);
    }

    return weights;
}

void run_time_step(std::vector<Vector> &points, std::vector<Vector> &velocity,
                   double *weights, std::default_random_engine &gen,
                   const double mass, const double eps, const double dt,
                   int frame_iter) {

    weights = optimal_transport_fluids(points, gen, weights);
    std::vector<Polygon> polygons = power_diagram_circle(points, weights);
    Vector g(0, -9.81, 0);

    for (size_t i = 0; i < polygons.size(); ++i) {
        Vector F_i_spring =
            1 / (eps * eps) * (compute_centroid(polygons[i]) - points[i]);

        Vector F_i = F_i_spring + mass * g;

        velocity[i] = velocity[i] + dt / mass * F_i;

        points[i] = points[i] + dt * velocity[i];
        if (points[i][1] < 0) {
            points[i][1] = 0;
        }
        if (points[i][1] > 1) {
            points[i][1] = 1;
        }

        if (points[i][0] < 0) {
            points[i][0] = 0;
        }
        if (points[i][0] > 1) {
            points[i][0] = 1;
        }
    }

    save_frame(polygons, "../final/out", frame_iter);
}

Polygon make_circle(const Vector &center, double radius) {
    constexpr int sides = 100;
    Polygon res;

    for (int i = 0; i < sides; ++i) {
        double x = radius * std::cos(2. * M_PI / sides * i);
        double y = radius * std::sin(2. * M_PI / sides * i);
        Vector point = Vector(x, y, 0) + center;
        point[0] = std::max(0., point[0]);
        point[0] = std::min(1., point[0]);
        point[1] = std::max(0., point[1]);
        point[1] = std::min(1., point[1]);
        res.add(point);
    }

    return res;
}

double compute_area_signed(const Polygon &polygon) {
    double area = 0.;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const Vector &cur_vertex = polygon[i];
        const Vector &prev_vertex =
            i > 0 ? polygon[i - 1] : polygon[polygon.size() - 1];

        area +=
            (prev_vertex[0] * cur_vertex[1] - cur_vertex[0] * prev_vertex[1]);
    }
    // Test compute_area with known polygons

    return area / 2.;
}

Vector compute_centroid(const Polygon &polygon) {

    double C_x = 0., C_y = 0.;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const Vector &cur_vertex = polygon[i];
        const Vector &prev_vertex =
            i > 0 ? polygon[i - 1] : polygon[polygon.size() - 1];

        C_x +=
            (cur_vertex[0] + prev_vertex[0]) *
            (prev_vertex[0] * cur_vertex[1] - cur_vertex[0] * prev_vertex[1]);

        C_y +=
            (cur_vertex[1] + prev_vertex[1]) *
            (prev_vertex[0] * cur_vertex[1] - cur_vertex[0] * prev_vertex[1]);
    }

    double area = compute_area_signed(polygon);

    // assert(area > 0);

    const double arr = 6 * area;

    return Vector(C_x / arr, C_y / arr, 0);
}

const Vector &Polygon ::operator[](size_t a) const { return vertices[a]; }
Vector &Polygon ::operator[](size_t a) { return vertices[a]; }

void Polygon::add(const Vector &__value) { vertices.push_back(__value); }

void Polygon::add(Vector &&__value) { vertices.push_back(__value); }

size_t Polygon::size() const { return vertices.size(); }

void save_svg(const std::vector<Polygon> &polygons, std::string filename,
              std::string fillcol) {
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" "
               "height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000),
                    (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

// Adds one frame of an animated svg file. frameid is the frame number (between
// 0 and nbframes-1). polygons is a list of polygons, describing the current
// frame. The polygon vertices are supposed to be in the range [0..1], and a
// canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons,
                       std::string filename, int frameid, int nbframes) {
    FILE *f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = "
                   "\"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000),
                    (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "	id = \"frame%u\"\n", frameid);
    fprintf(f, "	attributeName = \"display\"\n");
    fprintf(f, "	values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n	keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n	dur = \"5s\"\n");
    fprintf(f, "	begin = \"0s\"\n");
    fprintf(f, "	repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

int sgn(double val) { return (val > 0) - (val < 0); }

void save_frame(const std::vector<Polygon> &cells, std::string filename,
                int frameid) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W * H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W - 1., std::max(0., W * bminx));
        bminy = std::min(H - 1., std::max(0., H * bminy));
        bmaxx = std::max(W - 1., std::max(0., W * bmaxx));
        bmaxy = std::max(H - 1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 =
                        cells[i]
                            .vertices[(j + 1) % cells[i].vertices.size()][0] *
                        W;
                    double y1 =
                        cells[i]
                            .vertices[(j + 1) % cells[i].vertices.size()][1] *
                        H;
                    double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);

                    int sign = sgn(det);
                    if (prevSign == 0)
                        prevSign = sign;
                    else if (sign == 0)
                        sign = prevSign;
                    else if (sign != prevSign) {
                        isInside = false;
                        break;
                    }
                    prevSign = sign;
                    double edgeLen =
                        sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
                    double distEdge = std::abs(det) / edgeLen;
                    double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y0);
                    if (dotp < 0 || dotp > edgeLen * edgeLen)
                        distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    if (mindistEdge <= 2) {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 0;
                    } else {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 255;
                    }
                }
            }
        }
    }

    std::ostringstream os;

    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}
