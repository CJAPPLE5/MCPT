#ifndef HITTABLE_H
#define HITTABLE_H

#include "rtweekend.h"

#include "aabb.h"
#include "contants.h"

class material;
// #include "material.h"

struct hit_record
{
    point3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    double t = T_MAX;
    bool front_face;
    vec3 uv;
    int mesh_id;

    inline void set_face_normal(const ray &r, const vec3 &outward_normal)
    {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

// class hittable
// {
// public:
//     aabb aabb_box;
//     aabb get_box() { return aabb_box; }
//     virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const = 0;
//     virtual bool bounding_box(double time0, double time1, aabb &output_box) const = 0;

//     virtual double pdf_value(const vec3 &o, const vec3 &v) const
//     {
//         return 0.0;
//     }

//     virtual vec3 sample(const vec3 &o) const
//     {
//         return vec3(1, 0, 0);
//     }
// };

class hittable
{
public:
    hittable() {}

    hittable(const point3 &v0, const point3 &v1, const point3 &v2) : v0(v0), v1(v1), v2(v2), e1(v1 - v0), e2(v2 - v0), nf(normalize(cross(e1, e2)))
    {
        bbox.merge(v2);
        bbox.merge(v0);
        bbox.merge(v1);
        if (bbox.range(0) < EPSILON)
        {
            bbox.maximum.e[0] += EPSILON;
        }
        if (bbox.range(1) < EPSILON)
        {
            bbox.maximum.e[1] += EPSILON;
        }
        if (bbox.range(2) < EPSILON)
        {
            bbox.maximum.e[2] += EPSILON;
        }
    };

    hittable(const point3 &v0, const point3 &v1, const point3 &v2,
             shared_ptr<material> m) : v0(v0), v1(v1), v2(v2), e1(v1 - v0), e2(v2 - v0), nf(normalize(cross(e1, e2))), mat_ptr(m)
    {
        bbox.merge(v2);
        bbox.merge(v0);
        bbox.merge(v1);
        if (bbox.range(0) < EPSILON)
        {
            bbox.maximum.e[0] += EPSILON;
        }
        if (bbox.range(1) < EPSILON)
        {
            bbox.maximum.e[1] += EPSILON;
        }
        if (bbox.range(2) < EPSILON)
        {
            bbox.maximum.e[2] += EPSILON;
        }
    };

    hittable(const point3 &v0, const point3 &v1, const point3 &v2,
             const normal3 &n0, const normal3 &n1, const normal3 &n2,
             const vec3 &tex0, const vec3 &tex1, const vec3 &tex2,
             shared_ptr<material> m) : v0(v0), v1(v1), v2(v2), e1(v1 - v0), e2(v2 - v0),
                                       n0(n0), n1(n1), n2(n2),
                                       uv0(tex0), uv1(tex1), uv2(tex2),
                                       nf(normalize(cross(e1, e2))), mat_ptr(m)
    {
        bbox.merge(v0);
        bbox.merge(v1);
        bbox.merge(v2);
        if (bbox.range(0) < EPSILON)
        {
            bbox.maximum.e[0] += EPSILON;
        }
        if (bbox.range(1) < EPSILON)
        {
            bbox.maximum.e[1] += EPSILON;
        }
        if (bbox.range(2) < EPSILON)
        {
            bbox.maximum.e[2] += EPSILON;
        }
    }

    double getarea()
    {
        return abs(cross(e1, e2).length()) / 2;
    }
    bool
    hit(
        const ray &r, double t_min, double t_max, hit_record &rec) const;

    bool hit_(
        const ray &r, double t_min, double t_max, hit_record &rec) const;

    bool bounding_box(double time0, double time1, aabb &output_box) const;
    aabb get_bbox() { return bbox; }
    double pdf_value(const point3 &o, const vec3 &v) const;
    vec3 sample(const point3 &o) const;

    vec3 interpolatedNormal(double u, double v) const;
    vec3 interpolatedTexture(double u, double v) const;

public:
    point3 v0, v1, v2;
    normal3 nf;
    edge3 e1, e2;
    normal3 n0, n1, n2;
    shared_ptr<material> mat_ptr;
    vec3 uv0, uv1, uv2;
    aabb bbox;
    double tex_u, tex_v;
};

vec3 hittable::sample(const point3 &o) const
{
    double u = random_double(0, 1.0);
    double v = random_double(0, 1.0);
    return sqrt(v) * (1 - u) * v0 + u * sqrt(v) * v1 + (1 - sqrt(v)) * v2;
}

double hittable::pdf_value(const point3 &o, const vec3 &v) const
{
    return 1.0;
}

bool hittable::bounding_box(double time0, double time1, aabb &output_box) const
{
    output_box.merge(v0);
    output_box.merge(v1);
    output_box.merge(v2);
    return true;
}

vec3 hittable::interpolatedTexture(double u, double v) const
{
    return (1.0 - u - v) * uv0 + u * uv1 + v * uv2;
}

vec3 hittable::interpolatedNormal(double u, double v) const
{
    return normalize((1.0 - u - v) * n0 + u * n1 + v * n2);
}

bool hittable::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    vec3 P = cross(r.dir, e2);
    double determinant = dot(P, e1);
    if (determinant < EPSILON && determinant > -EPSILON)
    {
        return false;
    }

    double inv_determinant = 1.0 / determinant;

    vec3 T = r.orig - v0;
    double u = dot(P, T) * inv_determinant;
    if (u > 1.0 || u < 0)
    {
        return false;
    }

    vec3 Q = cross(T, e1);
    double v = dot(Q, r.dir) * inv_determinant;
    if (v > 1.0 || v < 0.0 || u + v > 1.0)
    {
        return false;
    }

    double t = dot(Q, e2) * inv_determinant;
    if (t <= t_min || t >= t_max)
    {
        return false;
    }

    rec.t = t;
    rec.p = r.at(rec.t);
    vec3 outward_normal = nf;
    if (!n0.near_zero())
    {
        outward_normal = interpolatedNormal(u, v);
    }
    // rec.set_face_normal(r, outward_normal);
    rec.normal = outward_normal;
    rec.mat_ptr = mat_ptr;
    rec.uv = interpolatedTexture(u, v);

    return true;
}

bool hittable::hit_(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    vec3 s = r.orig - v0;
    vec3 q = cross(r.dir, e2), rr = cross(s, e1);
    auto denominator = dot(q, e1);
    if (std::abs(denominator) < EPSILON)
    {
        return false;
    }
    vec3 w(dot(rr, e2), dot(q, s), dot(rr, r.dir));
    w /= denominator;
    if (w[0] < 1e-5)
    {
        return false;
    }
    if (w[0] < t_min || w[0] > t_max)
    {
        return false;
    }
    rec.t = w[0];
    w[0] = 1 - w[1] - w[2];

    rec.p = v0 * w[0] + v1 * w[1] + v2 * w[2];
    // std::cout << ret.point - old << std::endl;
    for (int i = 0; i < 3; ++i)
    {
        if (w[i] < -1e-5)
        {
            return false;
        }
    }
    rec.mat_ptr = mat_ptr;
    // rec.normal = n0 * w[0] + n1 * w[1] + n2 * w[2];
    rec.set_face_normal(r, n0 * w[0] + n1 * w[1] + n2 * w[2]);
    // rec.uv = uv[0] * w[0] + uv[1] * w[1] + uv[2] * w[2];
    // ret.obj = this;
    return true;
}

#endif
