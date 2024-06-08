#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <stack>
#include <stdio.h>
#include <string>
#include <unordered_set>
#include <vector>

#ifndef RAYTRACER_FILE
#define RAYTRACER_FILE 1
#include "raytracer.h"
#endif

#define INTERPOLATION 1 

#include "mesh.h"

void TriangleMesh::readOBJ(const char *obj) {

    char matfile[255];
    char grp[255];

    FILE *f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f))
            break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1],
                       &vec[2], &col[0], &col[1], &col[2]) == 6) {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);

            } else {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char *consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0,
                        &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0)
                    t.vtxi = vertices.size() + i0;
                else
                    t.vtxi = i0 - 1;
                if (i1 < 0)
                    t.vtxj = vertices.size() + i1;
                else
                    t.vtxj = i1 - 1;
                if (i2 < 0)
                    t.vtxk = vertices.size() + i2;
                else
                    t.vtxk = i2 - 1;
                if (j0 < 0)
                    t.uvi = uvs.size() + j0;
                else
                    t.uvi = j0 - 1;
                if (j1 < 0)
                    t.uvj = uvs.size() + j1;
                else
                    t.uvj = j1 - 1;
                if (j2 < 0)
                    t.uvk = uvs.size() + j2;
                else
                    t.uvk = j2 - 1;
                if (k0 < 0)
                    t.ni = normals.size() + k0;
                else
                    t.ni = k0 - 1;
                if (k1 < 0)
                    t.nj = normals.size() + k1;
                else
                    t.nj = k1 - 1;
                if (k2 < 0)
                    t.nk = normals.size() + k2;
                else
                    t.nk = k2 - 1;
                indices.push_back(t);
            } else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1,
                            &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2,
                                &offset);
                    if (nn == 3) {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0,
                                    &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (k0 < 0)
                            t.ni = normals.size() + k0;
                        else
                            t.ni = k0 - 1;
                        if (k1 < 0)
                            t.nj = normals.size() + k1;
                        else
                            t.nj = k1 - 1;
                        if (k2 < 0)
                            t.nk = normals.size() + k2;
                        else
                            t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true) {
                if (consumedline[0] == '\n')
                    break;
                if (consumedline[0] == '\0')
                    break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3) {
                    if (i0 < 0)
                        t2.vtxi = vertices.size() + i0;
                    else
                        t2.vtxi = i0 - 1;
                    if (i2 < 0)
                        t2.vtxj = vertices.size() + i2;
                    else
                        t2.vtxj = i2 - 1;
                    if (i3 < 0)
                        t2.vtxk = vertices.size() + i3;
                    else
                        t2.vtxk = i3 - 1;
                    if (j0 < 0)
                        t2.uvi = uvs.size() + j0;
                    else
                        t2.uvi = j0 - 1;
                    if (j2 < 0)
                        t2.uvj = uvs.size() + j2;
                    else
                        t2.uvj = j2 - 1;
                    if (j3 < 0)
                        t2.uvk = uvs.size() + j3;
                    else
                        t2.uvk = j3 - 1;
                    if (k0 < 0)
                        t2.ni = normals.size() + k0;
                    else
                        t2.ni = k0 - 1;
                    if (k2 < 0)
                        t2.nj = normals.size() + k2;
                    else
                        t2.nj = k2 - 1;
                    if (k3 < 0)
                        t2.nk = normals.size() + k3;
                    else
                        t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                } else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    } else {
                        nn =
                            sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (k0 < 0)
                                t2.ni = normals.size() + k0;
                            else
                                t2.ni = k0 - 1;
                            if (k2 < 0)
                                t2.nj = normals.size() + k2;
                            else
                                t2.nj = k2 - 1;
                            if (k3 < 0)
                                t2.nk = normals.size() + k3;
                            else
                                t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            } else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(f);
}

Intersect TriangleMesh::__intersect(const Ray &ray, int left, int right) const {
    Intersect current_intersect(false);
    for (int i = left; i < right; ++i) {
        const TriangleIndices &current_indicies = indices[i];
        const Vector &A = vertices[current_indicies.vtxi],
                     &normalA = normals[current_indicies.ni];
        const Vector &B = vertices[current_indicies.vtxj],
                     &normalB = normals[current_indicies.nj];
        const Vector &C = vertices[current_indicies.vtxk],
                     &normalC = normals[current_indicies.nk];
        const Vector e1 = B - A, e2 = C - A;

        // unnormalized normal
        Vector N = cross(e1, e2);
        double u_x_N = dot(ray.u, N);
        double beta = dot(e2, cross(A - ray.O, ray.u)) / u_x_N;
        double gamma = -dot(e1, cross(A - ray.O, ray.u)) / u_x_N;
        double alpha = 1 - beta - gamma;
        if (alpha <= 0. || alpha >= 1. || beta <= 0. || beta >= 1. || gamma <= 0. ||
            gamma >= 1.)
            continue;
        double t = dot(A - ray.O, N) / u_x_N;

        if(t > current_intersect.t) continue;

        Vector P = ray.O + t * ray.u;
        // Intersection normal (common behaviour normal)
        Vector Normal = alpha*normalA + beta*normalB + gamma* normalC;
        N.normalize();
        Normal.normalize();

        if (t < 0.)
            continue;

        if (t < current_intersect.t) {
            current_intersect.intersect = true;
            current_intersect.P = P;
#if INTERPOLATION == 1
            current_intersect.N = Normal;
#else
            current_intersect.N = N;
#endif

            current_intersect.t = t;
            current_intersect.object = this;
        }
    }

    return current_intersect;
}

Intersect TriangleMesh::intersect(const Ray &ray) const {
    int size = indices.size();

    double cur_dist, best_dist = std::numeric_limits<double>::max();

    if (!root->bbox->intersect(ray, cur_dist)) {
        return Intersect(false);
    }

    std::list<BVH *> my_stack;
    my_stack.push_back(root);

    Intersect best_intersect(false);
    Intersect curr_intersect;

    std::unordered_set<BVH *> visited;

    while (!my_stack.empty()) {
        auto bvh = my_stack.back();
        my_stack.pop_back();

        if (bvh->left_child != nullptr && bvh->right_child != nullptr) {

            if (bvh->left_child->bbox->intersect(ray, cur_dist)) {
                if (cur_dist < best_dist) {
                    best_dist = cur_dist;
                }
                my_stack.push_back(bvh->left_child);
            }

            if (bvh->right_child->bbox->intersect(ray, cur_dist)) {
                if (cur_dist < best_dist) {
                    best_dist = cur_dist;
                }
                my_stack.push_back(bvh->right_child);
            }
        } else {
            curr_intersect = __intersect(ray, bvh->left, bvh->right);
            if (curr_intersect.t < best_intersect.t) {
                best_intersect = curr_intersect;
            }
        }
    }
    return best_intersect;
}

void TriangleMesh::transform(double scale, const Vector &direction) {

    for (auto &vertex : vertices) {
        vertex = vertex * scale + direction;
    }

}

void BoundingBox::build_box(const TriangleMesh &triangle_mesh) {
    b_min = Vector(1e100, 1e100, 1e100);
    b_max = Vector(-1e100, -1e100, -1e100);

    for (const auto &vertex : triangle_mesh.vertices) {
        b_min[0] = std::min(b_min[0], vertex[0]);
        b_min[1] = std::min(b_min[1], vertex[1]);
        b_min[2] = std::min(b_min[2], vertex[2]);

        b_max[0] = std::max(b_max[0], vertex[0]);
        b_max[1] = std::max(b_max[1], vertex[1]);
        b_max[2] = std::max(b_max[2], vertex[2]);
    }
}

void BoundingBox::build_box(TriangleMesh &triangle_mesh, int left, int right) {
    b_min = Vector(1e128, 1e128, 1e128);
    b_max = Vector(-1e128, -1e128, -1e128);
    for (int i = left; i < right; ++i) {
        int idxA = triangle_mesh.indices[i].vtxi;
        int idxB = triangle_mesh.indices[i].vtxj;
        int idxC = triangle_mesh.indices[i].vtxk;

        Vector &vertex = triangle_mesh.vertices[idxA];

        b_min[0] = std::min(b_min[0], vertex[0]);
        b_min[1] = std::min(b_min[1], vertex[1]);
        b_min[2] = std::min(b_min[2], vertex[2]);

        b_max[0] = std::max(b_max[0], vertex[0]);
        b_max[1] = std::max(b_max[1], vertex[1]);
        b_max[2] = std::max(b_max[2], vertex[2]);

        Vector &vertex2 = triangle_mesh.vertices[idxB];

        b_min[0] = std::min(b_min[0], vertex2[0]);
        b_min[1] = std::min(b_min[1], vertex2[1]);
        b_min[2] = std::min(b_min[2], vertex2[2]);

        b_max[0] = std::max(b_max[0], vertex2[0]);
        b_max[1] = std::max(b_max[1], vertex2[1]);
        b_max[2] = std::max(b_max[2], vertex2[2]);

        Vector &vertex3 = triangle_mesh.vertices[idxC];
        b_min[0] = std::min(b_min[0], vertex3[0]);
        b_min[1] = std::min(b_min[1], vertex3[1]);
        b_min[2] = std::min(b_min[2], vertex3[2]);

        b_max[0] = std::max(b_max[0], vertex3[0]);
        b_max[1] = std::max(b_max[1], vertex3[1]);
        b_max[2] = std::max(b_max[2], vertex3[2]);
    }
}

bool BoundingBox::intersect(const Ray &ray, double &dist) const {

    const BoundingBox &bbox = *this;
    double tmin_x, tmax_x;
    double tmin_y, tmax_y;
    double tmin_z, tmax_z;
    tmin_x = (bbox.b_min[0] - ray.O[0]) / ray.u[0];
    tmax_x = (bbox.b_max[0] - ray.O[0]) / ray.u[0];
    tmin_y = (bbox.b_min[1] - ray.O[1]) / ray.u[1];
    tmax_y = (bbox.b_max[1] - ray.O[1]) / ray.u[1];
    tmin_z = (bbox.b_min[2] - ray.O[2]) / ray.u[2];
    tmax_z = (bbox.b_max[2] - ray.O[2]) / ray.u[2];

    if (tmin_x > tmax_x)
        std::swap(tmin_x, tmax_x);
    if (tmin_y > tmax_y)
        std::swap(tmin_y, tmax_y);
    if (tmin_z > tmax_z)
        std::swap(tmin_z, tmax_z);

    double tmin = std::max({tmin_x, tmin_y, tmin_z});
    double tmax = std::min({tmax_x, tmax_y, tmax_z});

    if (tmin <= tmax && tmax >= 0) {
        dist = tmin >= 0 ? tmin : tmax;
        return true;
    }

    return false;
}

Vector BVH::__compute_barycenter(TriangleMesh &triangle_mesh, int idx) {
    TriangleIndices &indices_coord = triangle_mesh.indices[idx];

    Vector &A = triangle_mesh.vertices[indices_coord.vtxi];
    Vector &B = triangle_mesh.vertices[indices_coord.vtxj];
    Vector &C = triangle_mesh.vertices[indices_coord.vtxk];

    return (A + B + C) / 3;
}
void BVH::__build_bvh(BVH *cur_bvh, TriangleMesh &triangle_mesh, int left_val,
                      int right_val) {

    cur_bvh->bbox = new BoundingBox(triangle_mesh, left_val, right_val);
    cur_bvh->left = left_val;
    cur_bvh->right = right_val;

    Vector diag = cur_bvh->bbox->b_max - cur_bvh->bbox->b_min;
    Vector middle_diag = cur_bvh->bbox->b_min + diag * 0.5;

    int long_axis = 0;
    if (diag[1] > diag[0]) {
        long_axis = 1;
    }
    if (diag[2] > diag[long_axis]) {
        long_axis = 2;
    }

    int pivot_index = left_val;

    for (int i = left_val; i < right_val; ++i) {
        Vector barycenter = __compute_barycenter(triangle_mesh, i);
        if (barycenter[long_axis] < middle_diag[long_axis]) {
            std::swap(triangle_mesh.indices[i],
                      triangle_mesh.indices[pivot_index]);
            ++pivot_index;
        }
    }

    if (pivot_index <= left_val || pivot_index >= (right_val - 1)) {
        return;
    }

    if ((right_val - left_val) <= 5) {
        return;
    }

    cur_bvh->left_child = new BVH();
    cur_bvh->right_child = new BVH();

    __build_bvh(cur_bvh->left_child, triangle_mesh, left_val, pivot_index);

    __build_bvh(cur_bvh->right_child, triangle_mesh, pivot_index, right_val);
}

void BVH::build_bvh(TriangleMesh &triangle_mesh) {
    int size = triangle_mesh.indices.size();

    __build_bvh(this, triangle_mesh, 0, size);
}

BVH::BVH(TriangleMesh &triangle_mesh) { build_bvh(triangle_mesh); }


