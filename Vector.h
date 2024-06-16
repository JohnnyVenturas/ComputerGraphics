#ifndef VECTOR_H
#define VECTOR_H
#include <random>
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
    Vector sample_direction(std::default_random_engine &engine,
                            std::uniform_real_distribution<double> &dist);

    // void operator=(const Vector &other) const;
    bool inside(const Vector &a, const Vector &b) const;

    double data[3];

    // static Vector intersect(const Vector &seg1, const Vector &seg2,
    // std::pair<Vector, Vector> &edge,
    // const Vector &pivot);
    //
    // static double compute_angle(const Vector &P, const Vector &A,
    // const Vector &B);
    //
    // static bool is_outwards(const Vector &N, const Vector &P1, const Vector
    // &P2, const Vector &P3);
    //
    bool inside(const Vector &bis1, const Vector &bis2);
};

Vector intersect(const Vector &seg1, const Vector &seg2, const Vector &bis1,
                 const Vector &bis2);

double cross_sign(const Vector &P, const Vector &A, const Vector &B);
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
std::ostream &operator<<(std::ostream &os, const Vector &v);

#endif
