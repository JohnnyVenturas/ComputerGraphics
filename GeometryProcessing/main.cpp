#include "Vector.h"
#include "mesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

double distance(const Vector &v1, const Vector &v2) {
    return std::sqrt(std::pow(v1[0] - v2[0], 2) + std::pow(v1[1] - v2[1], 2));
}

std::vector<int> getBoundaryIndices(const TriangleMesh &mesh) {
    std::vector<int> boundaryIndices;
    std::vector<int> edgeCount(mesh.vertices.size(), 0);

    for (const auto &triangle : mesh.indices) {
        edgeCount[triangle.vtxi]++;
        edgeCount[triangle.vtxj]++;
        edgeCount[triangle.vtxk]++;
    }

    for (int i = 0; i < edgeCount.size(); ++i) {
        if (edgeCount[i] == 1) {
            boundaryIndices.push_back(i);
        }
    }

    return boundaryIndices;
}

void initializeBoundaryVertices(TriangleMesh &mesh,
                                const std::vector<int> &boundaryIndices) {
    int n = boundaryIndices.size();
    double boundaryLength = 0.0;

    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        boundaryLength += distance(mesh.vertices[boundaryIndices[i]],
                                   mesh.vertices[boundaryIndices[next]]);
    }

    double cs = 0.0;
    for (int i = 0; i < n; ++i) {
        double theta = 2 * M_PI * cs / boundaryLength;
        mesh.vertices[boundaryIndices[i]][0] = std::cos(theta);
        mesh.vertices[boundaryIndices[i]][1] = std::sin(theta);
        mesh.vertices[boundaryIndices[i]][2] = 0.0;

        int next = (i + 1) % n;
        cs += distance(mesh.vertices[boundaryIndices[i]],
                       mesh.vertices[boundaryIndices[next]]);
    }
}

void updateInternalVertices(TriangleMesh &mesh,
                            const std::vector<int> &boundaryIndices,
                            int nbIter) {
    int numVertices = mesh.vertices.size();
    for (int iter = 0; iter < nbIter; ++iter) {
        std::vector<Vector> newVertices = mesh.vertices;

        for (int i = 0; i < numVertices; ++i) {
            if (std::find(boundaryIndices.begin(), boundaryIndices.end(), i) ==
                boundaryIndices.end()) {
                Vector newPos(0.0, 0.0, 0.0);
                int K = 0;
                for (const auto &triangle : mesh.indices) {
                    if (triangle.vtxi == i || triangle.vtxj == i ||
                        triangle.vtxk == i) {
                        if (triangle.vtxi != i)
                            newPos += mesh.vertices[triangle.vtxi];
                        if (triangle.vtxj != i)
                            newPos += mesh.vertices[triangle.vtxj];
                        if (triangle.vtxk != i)
                            newPos += mesh.vertices[triangle.vtxk];
                        K += 2;
                    }
                }
                newPos = newPos / K;
                newVertices[i] = newPos;
            }
        }

        mesh.vertices = newVertices;
    }
}

int main() { 
    TriangleMesh mesh("goethe.obj", ObjectProperties()); 
    std::vector<int> boundaryIndices = getBoundaryIndices(mesh);

    initializeBoundaryVertices(mesh, boundaryIndices);

    int nbIter = 100;
    updateInternalVertices(mesh, boundaryIndices, nbIter);

    return 0;
}

