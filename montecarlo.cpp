#include <cmath>
#include <iostream>
#include <random>

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0,1);

void box_muller(double stdev, double &x, double &y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
    y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

double F(double x, double y, double z) {
  return cos(x * y * z);
}

double p(double x, double y, double z, double stdev) {
  const double constant = 1 / (stdev) * 1 / (sqrt(2 * M_PI));
  const double constant3 = constant * constant * constant;
  return constant3 * exp(-(x*x + y *y + z*z) / (2 * stdev * stdev));
}

double compute_integral(double (*F)(double, double,double), double (*p)(double, double, double ,double), double stdev, int n) {
  double sum = 0.;
  double x, y,z, tmp;
  for(int i = 0; i < n; ++i) {
    box_muller(stdev, x,y);
    box_muller(stdev, z,tmp);
    if(x > M_PI / 2 || x < -M_PI / 2) continue;
    if(y > M_PI / 2 || y < -M_PI / 2) continue;
    if(z > M_PI / 2 || z < -M_PI / 2) continue;
    sum += F(x,y,z) / p(x,y,z,stdev) * 1./n;
  }

  return sum;

}


int main() {

  std::cout << compute_integral(F, p, 1, 1000000) << "\n";

}
