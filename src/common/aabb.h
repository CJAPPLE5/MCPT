#ifndef AABB_H
#define AABB_H

#include "rtweekend.h"

class aabb
{
public:
    aabb()
    {
        maximum = vec3(-T_MAX, -T_MAX, -T_MAX);
        minimum = vec3(T_MAX, T_MAX, T_MAX);
    }
    aabb(const point3 &a, const point3 &b)
    {
        minimum = min_vector(minimum, a);
        minimum = min_vector(minimum, b);
        maximum = max_vector(maximum, a);
        maximum = max_vector(maximum, b);
    }

    point3 min() const { return minimum; }
    point3 max() const { return maximum; }

    double range(int axis)
    {
        if (axis == 0)
            return maximum.x() - minimum.x();
        else if (axis == 1)
            return maximum.y() - minimum.y();
        else
            return maximum.z() - minimum.z();
    }

    void merge(const point3 &p)
    {
        minimum = min_vector(minimum, p);
        maximum = max_vector(maximum, p);
    }

    point3 get_center()
    {
        return (maximum + minimum) / 2;
    }

    void combine(const aabb &box_)
    {
        merge(box_.min());
        merge(box_.max());
    }

    bool hit(const ray &r, double t_min, double t_max) const
    {
        for (int a = 0; a < 3; a++)
        {
            if (r.direction()[a] == 0)
            {
                continue;
            }
            auto t0 = fmin((minimum[a] - r.origin()[a]) / (r.direction()[a]),
                           (maximum[a] - r.origin()[a]) / (r.direction()[a]));
            auto t1 = fmax((minimum[a] - r.origin()[a]) / (r.direction()[a]),
                           (maximum[a] - r.origin()[a]) / (r.direction()[a]));
            t_min = fmax(t0, t_min);
            t_max = fmin(t1, t_max);
            if (t_max <= t_min)
                return false;
        }
        return true;

        // vec3 o = r.orig, d = r.dir;
        // float t0 = T_MIN, t1 = T_MAX;

        // if (std::fabs(d.x()) > EPSILON)
        // {
        //     t0 = fmax(t0, ((d.x() > 0 ? minimum.x() : maximum.x()) - o.x()) / d.x());
        //     t1 = fmin(t1, ((d.x() > 0 ? maximum.x() : minimum.x()) - o.x()) / d.x());
        // }
        // else if (o.x() < minimum.x() || o.x() > maximum.x())
        //     return false;

        // if (std::fabs(d.y()) > EPSILON)
        // {
        //     t0 = fmax(t0, ((d.y() > 0 ? minimum.y() : maximum.y()) - o.y()) / d.y());
        //     t1 = fmin(t1, ((d.y() > 0 ? maximum.y() : minimum.y()) - o.y()) / d.y());
        // }
        // else if (o.y() < minimum.y() || o.y() > maximum.y())
        //     return false;

        // if (std::fabs(d.z()) > EPSILON)
        // {
        //     t0 = fmax(t0, ((d.z() > 0 ? minimum.z() : maximum.z()) - o.z()) / d.z());
        //     t1 = fmin(t1, ((d.z() > 0 ? maximum.z() : minimum.z()) - o.z()) / d.z());
        // }
        // else if (o.z() < minimum.z() || o.z() > maximum.z())
        //     return false;

        // return t0 <= t1;

        // vec3 t1 = (minimum - r.orig) / r.dir;
        // vec3 t2 = (maximum - r.orig) / r.dir;
        // vec3 max_t = max_vector(t1, t2);
        // vec3 min_t = min_vector(t1, t2);
        // float total_mint = std::max(min_t[0], std::max(min_t[1], min_t[2]));
        // float total_maxt = std::min(max_t[0], std::min(max_t[1], max_t[2]));
        // return total_maxt >= 0 && total_mint <= total_maxt;
    }

    double area() const
    {
        auto a = maximum.x() - minimum.x();
        auto b = maximum.y() - minimum.y();
        auto c = maximum.z() - minimum.z();
        return 2 * (a * b + b * c + c * a);
    }

    int longest_axis() const
    {
        auto a = maximum.x() - minimum.x();
        auto b = maximum.y() - minimum.y();
        auto c = maximum.z() - minimum.z();
        if (a > b && a > c)
            return 0;
        else if (b > c)
            return 1;
        else
            return 2;
    }

public:
    point3 minimum;
    point3 maximum;
};

aabb surrounding_box(aabb box0, aabb box1)
{
    vec3 small(fmin(box0.min().x(), box1.min().x()),
               fmin(box0.min().y(), box1.min().y()),
               fmin(box0.min().z(), box1.min().z()));

    vec3 big(fmax(box0.max().x(), box1.max().x()),
             fmax(box0.max().y(), box1.max().y()),
             fmax(box0.max().z(), box1.max().z()));

    return aabb(small, big);
}

#endif
