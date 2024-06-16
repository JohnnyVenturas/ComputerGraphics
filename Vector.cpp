#include "Vector.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <random>

double epsilon = 1e-5;
Vector::Vector(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

std::ostream &operator<<(std::ostream &os, const Vector &v) {

    os << "(";
    os << std::fixed << std::setw(10) << std::setprecision(6)
       << std::setfill(' ') << v[0];
    os << ", ";
    os << std::fixed << std::setw(10) << std::setprecision(6)
       << std::setfill(' ') << v[1];
    os << ", ";

    os << std::fixed << std::setw(10) << std::setprecision(6)
       << std::setfill(' ') << v[2];
    os << ")\n";

    return os;
}

// Vector Vector ::intersect(const Vector &seg1, const Vector &seg2,
// std::pair<Vector, Vector> &edge,
// const Vector &pivot) {
//
//// ln stands for line
//// seg stand for segment
//// ln1, ln2 corresponds to the pair(ln1, ln2) meaning from ln1 to ln2
//// same for seg1, seg2
//
// const Vector &ln1 = edge.first;
// const Vector &ln2 = edge.second;
//
// Vector N = Vector(ln2[1] - ln1[1], ln2[0] - ln1[0], ln2[2] - ln1[2]);
//
// if (!Vector::is_outwards(N, ln1, ln2, pivot)) {
// N = -N;
//}
//
//// make sure N is outwards
//
// double t = dot(ln1 - seg1, N) / dot(seg2 - seg1, N);
//
// Vector P = seg1 + t * (seg2 - seg1);
//
// if (t < 0 || t > 1.)
// throw "This shall not happen.";
//
// return P;
//}
//
// bool Vector::is_outwards(const Vector &N, const Vector &P1, const Vector &P2,
// const Vector &P3) {
//
// double angle1 = dot(P2 - P1, N);
// double angle2 = dot(P3 - P2, N);
//
// if (angle1 == angle2) {
// return false;
//}
//
// return true;
//}
//
// bool Vector ::inside(const std::pair<Vector, Vector> &edge,
// const Vector &pivot) {
//
// const Vector &P = *this;
// const Vector &u = edge.first, &v = edge.second;
//
// Vector N = Vector(v[1] - u[1], v[0] - u[0], v[2] - u[2]);
// if (!Vector::is_outwards(N, u, v, pivot)) {
// N = -N;
//}
//
// if (dot(P - u, N) < 0.)
// return true;
//
// return false;
//}

bool Vector ::inside(const Vector &bis1, const Vector &bis2) const {
    const Vector &p = *this;

    return (p - bis1).norm2() <= (p - bis2).norm2();
}

Vector intersect(const Vector &seg1, const Vector &seg2, const Vector &p0,
                 const Vector &p1) {

    const Vector M = (p0 + p1) / 2;

    double t = dot(M - seg1, p1 - p0) / dot(seg2 - seg1, p1 - p0);

    return seg1 + t * (seg2 - seg1);
}

Vector Vector ::sample_direction(std::default_random_engine &gen,
                                 std::uniform_real_distribution<double> &dist) {

    Vector normal = *this;
    double r1 = dist(gen), r2 = dist(gen);

    Vector T1 = normal;

    int pos = 0;

    if (std::abs(T1.data[0]) <= std::abs(T1.data[1]) &&
        std::abs(T1.data[0]) <= std::abs(T1.data[2])) {
        pos = 0;
    } else if (std::abs(T1.data[1]) <= std::abs(T1.data[2]) &&
               std::abs(T1.data[1]) <= std::abs(T1.data[0])) {
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
    return std::abs(data[0] - other.data[0]) < epsilon &&
           std::abs(data[1] - other.data[1]) < epsilon &&
           std::abs(data[2] - other.data[2]) < epsilon;
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

