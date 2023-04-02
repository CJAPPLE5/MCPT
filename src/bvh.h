#ifndef BVH_H
#define BVH_H

#include "rtweekend.h"

#include "hittable.h"

#include <algorithm>

class bvh_node
{
public:
    bvh_node();

    bvh_node(const std::vector<shared_ptr<hittable>> &src_objects);

    bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec) const;

    bool bounding_box(double time0, double time1, aabb &output_box) const;
    std::vector<shared_ptr<hittable>> triangles;

    shared_ptr<bvh_node> left;
    shared_ptr<bvh_node> right;
    aabb box;
};

bool compareX(const shared_ptr<hittable> &h0, const shared_ptr<hittable> &h1)
{
    return h0->get_bbox().get_center().x() < h1->get_bbox().get_center().x();
}

bool compareY(const shared_ptr<hittable> &h0, const shared_ptr<hittable> &h1)
{
    return h0->get_bbox().get_center().y() < h1->get_bbox().get_center().y();
}

bool compareZ(const shared_ptr<hittable> &h0, const shared_ptr<hittable> &h1)
{
    return h0->get_bbox().get_center().z() < h1->get_bbox().get_center().z();
}

inline bool box_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis)
{
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    return box_a.min().e[axis] < box_b.min().e[axis];
}

bvh_node::bvh_node(
    const std::vector<shared_ptr<hittable>> &src_objects)
{
    for (int i = 0; i < src_objects.size(); i++)
    {
        // aabb temp_box;
        box.combine(src_objects[i]->get_bbox());
    }
    if (src_objects.size() <= BVH_MAX)
    {
        this->triangles = src_objects;
        left = right = nullptr;
    }
    else
    {
        double x = box.range(0), y = box.range(1), z = box.range(2);
        std::vector<shared_ptr<hittable>> tempTriangles = src_objects;
        int axis = random_int(0, 2);

        if (x >= y && x >= z)
            std::sort(tempTriangles.begin(), tempTriangles.end(), compareX);
        else if (y >= x && y >= z)
            std::sort(tempTriangles.begin(), tempTriangles.end(), compareY);
        else
            std::sort(tempTriangles.begin(), tempTriangles.end(), compareZ);
        auto middle = tempTriangles.begin() + (tempTriangles.size() / 2);
        std::vector<shared_ptr<hittable>> leftTriangles(tempTriangles.begin(), middle);
        std::vector<shared_ptr<hittable>> rightTriangles(middle, tempTriangles.end());
        left = make_shared<bvh_node>(leftTriangles);
        right = make_shared<bvh_node>(rightTriangles);
    }
}

bool bvh_node::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    if (!box.hit(r, t_min, t_max))
    {
        // std::cout << "reject" << std::endl;
        return false;
    }
    if (triangles.size() > 0)
    {
        bool hit_any = false;
        double t_sofar = t_max;
        for (int i = 0; i < triangles.size(); i++)
        {
            if (triangles[i]->hit(r, t_min, t_sofar, rec))
            {
                hit_any = true;
                t_sofar = rec.t;
            }
        }
        return hit_any;
    }
    bool hit_left = false, hit_right = false;
    if (left != nullptr)
        hit_left = left->hit(r, t_min, t_max, rec);
    if (right != nullptr)
        hit_right = right->hit(r, t_min, fmin(rec.t, t_max), rec);

    return hit_left || hit_right;
}

bool bvh_node::bounding_box(double time0, double time1, aabb &output_box) const
{
    output_box = box;
    return true;
}

#endif
