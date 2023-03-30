#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "hittable.h"

// class triangle : public hittable
// {
// public:
//     triangle() {}

//     triangle(const point3 &v0, const point3 &v1, const point3 &v2,
//              shared_ptr<material> m) : v0(v0), v1(v1), v2(v2), e1(v1 - v0), e2(v2 - v0), nf(normalize(cross(e1, e2))), mat_ptr(m), bbox(v0, v1) { bbox.merge(v2); };

//     triangle(const point3 &v0, const point3 &v1, const point3 &v2,
//              const normal3 &n1, const normal3 &n2, const normal3 &n3,
//              shared_ptr<material> m) : v0(v0), v1(v1), v2(v2), e1(v1 - v0), e2(v2 - v0), n0(n0), n1(n1), n2(n2),
//                                        nf(normalize(cross(e1, e2))), mat_ptr(m), bbox(v0, v1) { bbox.merge(v1); }

//     virtual bool hit(
//         const ray &r, double t_min, double t_max, hit_record &rec) const override;

//     virtual bool bounding_box(double time0, double time1, aabb &output_box) const override;
//     virtual double pdf_value(const point3 &o, const vec3 &v) const override;
//     // virtual vec3 random(const point3 &o) const override;

//     vec3 interpolatedNormal(double u, double v) const;

// public:
//     point3 v0, v1, v2;
//     normal3 nf;
//     edge3 e1, e2;
//     normal3 n0, n1, n2;
//     shared_ptr<material> mat_ptr;
//     aabb bbox;
// };

// double triangle::pdf_value(const point3 &o, const vec3 &v) const
// {
//     return 1.0;
// }

// bool triangle::bounding_box(double time0, double time1, aabb &output_box) const
// {
//     output_box.merge(v0);
//     output_box.merge(v1);
//     output_box.merge(v2);
//     return true;
// }

// vec3 triangle::interpolatedNormal(double u, double v) const
// {
//     return normalize((1.0 - u - v) * n0 + u * n1 + v * n2);
// }

// bool triangle::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
// {
//     vec3 P = cross(r.dir, e2);
//     double determinant = dot(P, e1);
//     if (determinant < EPSILON && determinant > -EPSILON)
//     {
//         return false;
//     }

//     double inv_determinant = 1.0 / determinant;

//     vec3 T = r.orig - v0;
//     double u = dot(P, T) * inv_determinant;
//     if (u > 1.0 || u < 0.0)
//     {
//         return false;
//     }

//     vec3 Q = cross(T, e1);
//     double v = dot(Q, r.dir) * inv_determinant;
//     if (v > 1.0 || v < 0.0 || u + v > 1.0)
//     {
//         return false;
//     }

//     double t = dot(Q, e2) * inv_determinant;
//     if (t <= 0.0)
//     {
//         return false;
//     }

//     rec.t = t;
//     rec.p = r.at(rec.t);
//     vec3 outward_normal = nf;
//     if (!n0.near_zero())
//     {
//         outward_normal = interpolatedNormal(u, v);
//     }
//     rec.set_face_normal(r, outward_normal);
//     rec.mat_ptr = mat_ptr;

//     return true;
// }

#endif