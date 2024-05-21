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
double dot(const Vector &a, const Vector &b);
Vector cross(const Vector &a, const Vector &b);
