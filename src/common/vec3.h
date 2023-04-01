#ifndef VEC3_H
#define VEC3_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include <cmath>
#include <iostream>

using std::fabs;
using std::sqrt;

class vec3
{
public:
    vec3() : e{0, 0, 0} {}
    vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double &operator[](int i) { return e[i]; }

    vec3 &operator+=(const vec3 &v)
    {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vec3 &operator*=(const double t)
    {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    vec3 &operator/=(const double t)
    {
        return *this *= 1 / t;
    }

    double length() const
    {
        return sqrt(length_squared());
    }

    double length_squared() const
    {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    bool near_zero() const
    {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }

    inline static vec3 random()
    {
        return vec3(random_double(), random_double(), random_double());
    }

    inline static vec3 random(double min, double max)
    {
        return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
    }

public:
    double e[3];
};

// Type aliases for vec3
using point3 = vec3; // 3D point
using color = vec3;  // RGB color
using edge3 = vec3;
using normal3 = vec3;

// vec3 Utility Functions

inline std::ostream &operator<<(std::ostream &out, const vec3 &v)
{
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 abs(const vec3 &u)
{
    return vec3(abs(u.e[0]), abs(u.e[1]), abs(u.e[2]));
}

inline vec3 operator+(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator/(const vec3 &u, const vec3 &v)
{

    return vec3((v.e[0] != 0) ? u.e[0] / v.e[0] : u.e[0],
                (v.e[1] != 0) ? u.e[1] / v.e[1] : u.e[1],
                (v.e[2] != 0) ? u.e[2] / v.e[2] : u.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3 &v)
{
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3 &v, double t)
{
    return t * v;
}

inline vec3 operator/(vec3 v, double t)
{
    return (1 / t) * v;
}

inline double dot(const vec3 &u, const vec3 &v)
{
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 min_vector(const vec3 &a, const vec3 &b)
{
    return vec3(fmin(a.x(), b.x()), fmin(a.y(), b.y()), fmin(a.z(), b.z()));
}

inline vec3 max_vector(const vec3 &a, const vec3 &b)
{
    return vec3(fmax(a.x(), b.x()), fmax(a.y(), b.y()), fmax(a.z(), b.z()));
}

inline vec3 unit_vector(vec3 v)
{
    return v / v.length();
}

inline vec3 random_in_unit_disk()
{
    while (true)
    {
        auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
        if (p.length_squared() >= 1)
            continue;
        return p;
    }
}

inline vec3 random_in_unit_sphere()
{
    while (true)
    {
        auto p = vec3::random(-1, 1);
        if (p.length_squared() >= 1)
            continue;
        return p;
    }
}

inline vec3 normalize(const vec3 &v)
{
    return v / v.length();
}

inline vec3 random_unit_vector()
{
    return unit_vector(random_in_unit_sphere());
}

inline vec3 random_in_hemisphere(const vec3 &normal)
{
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline void sampleHemisphere(const double exp, double &theta, double &phi)
{
    // theta = std::acos(std::pow(random_double(), 1.0f / (exp + 1.0f)));
    // phi = random_double() * PI * 2.0f;
    theta = std::pow(1 - random_double(), 1.0 / (exp + 1.0));
}

inline void calculateTangentSpace(const vec3 &normal, vec3 &x, vec3 &y)
{
    vec3 z = normalize(normal);
    vec3 a = (fabs(z.x()) > 0.9) ? vec3(0, 1, 0) : vec3(1, 0, 0);
    x = normalize(cross(z, a));
    y = normalize(cross(z, x));
}

inline vec3 random_hemisphere(const vec3 &wi, const vec3 &normal, double ns, double &pdf, const double &pd)
{
    // double theta = std::acos(std::pow(1 - random_double(), 1.0 / (ns)));
    float cos_theta = std::sqrt(1 - random_double());
    float sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    double phi = 2 * PI * random_double();
    vec3 x, y;
    calculateTangentSpace(normal, x, y);
    vec3 direction = cos_theta * normal + sin_theta * std::cos(phi) * x + sin_theta * std::sin(phi) * y;
    pdf = pd * cos_theta / PI;
    if (ns > 1)
    {
        vec3 half = normalize(wi + direction);
        double tmp_cos1 = dot(half, normal);
        double tmp_cos2 = dot(half, wi);
        pdf += (1 - pd) * (ns + 2) / TWO_PI * std::pow(tmp_cos1, ns + 1) / (4 * tmp_cos2);
    }

    return direction;
}

inline vec3 random_hemisphere_specular(const vec3 &wi, const vec3 &normal, double ns, double &pdf, double &pd)
{
    double cos_theta = std::pow(1 - random_double(), 1.0f / (ns + 2.0f));
    double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    double phi = TWO_PI * random_double();

    vec3 x, y;
    calculateTangentSpace(normal, x, y);
    vec3 local = cos_theta * normal + sin_theta * std::cos(phi) * x + sin_theta * std::sin(phi) * y;
    double temp_cos = dot(local, wi);
    vec3 wo = local * 2 * temp_cos - wi;
    pdf = (1 - pd) * (ns + 2) / TWO_PI * std::pow(cos_theta, ns + 1);
    pdf /= (4 * temp_cos);
    if (dot(normal, wo) <= 0)
    {
        pdf = 0;
    }
    else
    {
        pdf += pd * dot(normal, wo) / PI;
    }
    return wo;
}

inline vec3 reflect(const vec3 &v, const vec3 &n)
{
    return v - 2 * dot(v, n) * n;
}

inline vec3 refract(const vec3 &uv, const vec3 &n, double etai_over_etat)
{
    auto cos_theta = fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

class matrix
{
public:
    double e[4][4];

    matrix() {}

    matrix(const double **e_)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                e[i][j] = e_[i][j];
            }
        }
    }
    vec3 mlt(const vec3 &v) const
    {
        vec3 result;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                result[i] += e[i][j] * v.e[j];
            }
        }
        return result;
    }

    matrix transpose() const
    {
        matrix result;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                result.e[i][j] = e[j][i];
            }
        }
        return result;
    }
};

#endif
