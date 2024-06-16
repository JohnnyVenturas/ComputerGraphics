#include "Vector.h"
#include "polygon.h"
#include "raytracer.h"
#include <random>
#include <chrono>
#include <vector>

std::default_random_engine
    gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
int main() {
    int W = 512;
    int H = 512;

    // draw_image(W, H);

    draw_random_polygon(gen);

//    test_area();

    return 0;
}
