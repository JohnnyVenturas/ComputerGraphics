Here is the revised README:

---

# Computer Graphics Raytracer

## Author
Bogdan Budura

## Date
May 2024

## Overview
This project is a raytracer that simulates light behavior to create realistic images. It implements various features such as diffuse and mirror surfaces, indirect lighting, anti-aliasing, and complex object handling using triangle meshes and Bounding Volume Hierarchies (BVH).

## Features

### Diffuse and Mirror Surfaces
- **Scene Setup**: A scene with six colored walls and three spheres (one reflective and two refractive).
- **Implementation**: Light interactions are simulated using 6 light bounces and one sample per pixel, achieving an execution time of 35ms. Enhanced with Fresnel Law for partial reflection and refraction, averaging 1000 rays per pixel.

### Indirect Lighting
- **Implementation**: Simulates diffusive surface reflections by randomly sampling light directions and recursively computing light incidence over multiple bounces.
- **Anti-Aliasing**: Addresses visual aberrations by sampling intersection points within the pixel surface, averaging 1000 rays per pixel.

### Cat Mesh
- **Complex Objects**: Handles complex objects using triangle meshes for detailed approximation.
- **Bounding Volume Hierarchies (BVH)**: Optimizes performance by enclosing objects in hierarchical bounding boxes to reduce computationally expensive intersection tests.

## Results
- **Diffuse and Mirror Surfaces**: Initial rendering in 35ms, enhanced scene with Fresnel Law rendered in approximately 13 seconds.
- **Indirect Lighting**: Rendered with 1000 rays per pixel in 58.5 seconds; combined with anti-aliasing in 1 minute and 5 seconds.
- **Cat Mesh**: Initial rendering took over 40 minutes; optimized with BVH to 9.65 seconds with 64 samples per pixel.

## Conclusion
The raytracer project provided valuable lessons in programming and problem-solving, with significant challenges in debugging and optimizing performance. The implementation of BVH was particularly difficult but resulted in substantial performance improvements and high-quality renders.

---

