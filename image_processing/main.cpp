#define _CRT_SECURE_NO_WARNINGS 1
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <Vector.h> 

Vector random_direction(std::default_random_engine &engine,
                        std::uniform_real_distribution<double> &dist) {
    double r1 = dist(engine);
    double r2 = dist(engine);
    double x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double z = 1 - 2 * r2;
    return Vector(x, y, z);
}

void sliced_optimal_transport(std::vector<Vector> &input_image,
                              std::vector<Vector> &model_image, int n,
                              int nbiter) {
    std::default_random_engine engine;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int iter = 0; iter < nbiter; ++iter) {
        Vector v = random_direction(engine, dist);

        std::vector<std::pair<double, int>> projI(n), projM(n);
        for (int i = 0; i < n; ++i) {
            projI[i] = {dot(input_image[i], v), i};
            projM[i] = {dot(model_image[i], v), i};
        }

        std::sort(projI.begin(), projI.end());
        std::sort(projM.begin(), projM.end());

        for (int i = 0; i < n; ++i) {
            int idxI = projI[i].second;
            int idxM = projM[i].second;
            input_image[idxI] += (model_image[idxM] - input_image[idxI]) * v;
        }
    }
}

int main() {
    int W, H, C;

    unsigned char *image =
        stbi_load("8733654151_b9422bb2ec_k.jpg", &W, &H, &C, STBI_rgb);
    if (!image) {
        std::cerr << "Error loading input image." << std::endl;
        return -1;
    }

    unsigned char *model_image_data =
        stbi_load("model.jpg", &W, &H, &C, STBI_rgb);
    if (!model_image_data) {
        std::cerr << "Error loading model image." << std::endl;
        stbi_image_free(image);
        return -1;
    }

    std::vector<Vector> input_image(W * H);
    std::vector<Vector> model_image(W * H);
    for (int i = 0; i < W * H; ++i) {
        input_image[i] =
            Vector(image[3 * i], image[3 * i + 1], image[3 * i + 2]);
        model_image[i] =
            Vector(model_image_data[3 * i], model_image_data[3 * i + 1],
                   model_image_data[3 * i + 2]);
    }

    int nbiter = 100;
    sliced_optimal_transport(input_image, model_image, W * H, nbiter);

    std::vector<double> image_double(W * H * 3);
    for (int i = 0; i < W * H; ++i) {
        image_double[3 * i] = input_image[i][0];
        image_double[3 * i + 1] = input_image[i][1];
        image_double[3 * i + 2] = input_image[i][2];
    }

    std::vector<unsigned char> image_result(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            image_result[(i * W + j) * 3 + 0] =
                image_double[(i * W + j) * 3 + 0] * 0.5;
            image_result[(i * W + j) * 3 + 1] =
                image_double[(i * W + j) * 3 + 1] * 0.3;
            image_result[(i * W + j) * 3 + 2] =
                image_double[(i * W + j) * 3 + 2] * 0.2;
        }
    }
    stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

    return 0;
}

