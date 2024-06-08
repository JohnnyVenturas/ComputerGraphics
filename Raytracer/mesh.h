#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#ifndef RAYTRACER_FILE
#define RAYTRACER_FILE 1
#include "raytracer.h"
#endif

class TriangleIndices {
  public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1,
                    int nj = -1, int nk = -1, int uvi = -1, int uvj = -1,
                    int uvk = -1, int group = -1, bool added = false)
        : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk),
          ni(ni), nj(nj), nk(nk), group(group){};

    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class TriangleMesh;

class BoundingBox {
  public:
    
    Vector b_min, b_max;
    BoundingBox() {}

    BoundingBox(const Vector &b_min, const Vector &b_max)
        : b_min(b_min), b_max(b_max) {}

    BoundingBox(const TriangleMesh &triangle_mesh) {
        build_box(triangle_mesh);
    }

    BoundingBox(TriangleMesh &triangle_mesh, int left, int right) {
        build_box(triangle_mesh, left, right);
    }


    void build_box(const TriangleMesh &triangle_mesh);
    void build_box(TriangleMesh &triangle_mesh, int left, int right);
    
    bool intersect(const Ray &ray, double &dist) const;

};

class TriangleMesh;
class BVH {
private:
  Vector __compute_barycenter(TriangleMesh &triangle_mesh, int idx);
  void __build_bvh(BVH *cur_bvh, TriangleMesh &triangle_mesh, int left_val, int right_val);

public:
  BoundingBox * bbox;
  BVH *left_child = nullptr, *right_child = nullptr;
  int left, right;

  BVH() {} 

  BVH(TriangleMesh &triangle_mesh);

  void build_bvh(TriangleMesh &triangle_mesh);
  bool intersect (const Ray &ray);

};

class TriangleMesh : public virtual Geometry {
private:
    Intersect __intersect(const Ray &ray,int left, int right) const;

public:
    ~TriangleMesh() {}
    TriangleMesh(){};

    TriangleMesh(const char *obj, const ObjectProperties &object_properties)
        : Geometry(object_properties) {
        readOBJ(obj);

    }

    BVH *root;

    void readOBJ(const char *obj);

    void transform(double scale, const Vector &direction);

    virtual Intersect intersect(const Ray &ray) const;

    bool bbox_intersect(const Ray &ray, const BoundingBox &bbox) const;
    bool bvh_box_intersect(const Ray &ray) ;

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
};

