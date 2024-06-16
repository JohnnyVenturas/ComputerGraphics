#include <cmath>
#include <random>
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include <iostream>

#include "Scene.h"
#include "mesh.h"
#include "raytracer.h"

// static std::random_device rd;

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
    TriangleMesh *cat_mesh = new TriangleMesh("../cat/cat.obj", reflective);

    cat_mesh->transform(0.6, Vector(0, -10, -10));
    cat_mesh->root = new BVH(*cat_mesh);
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
                Ray camera_ray(camera, x, y, ANTIALIASING, gen, dist);

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

    // std::cout <<"Sanity check\n";
    // std::cout << static_cast<Geometry *>(cat_mesh) << " " << cat_mesh ->
    // object_properties.material_type<< "\n"; std::cout <<"End\n";

    stbi_write_png("image.png", W, H, 3, &image[0], 0);
}

