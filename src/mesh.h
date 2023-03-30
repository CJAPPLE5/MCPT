#ifndef MESH_H
#define MESH_H

#include "hittable.h"
#include "bvh.h"
class Mesh
{
public:
    double area;
    shared_ptr<bvh_node> mesh_bvh;
    // std::vector<shared_ptr<triangle>> triangles;
    std::vector<double> tri_areas;
    std::vector<shared_ptr<hittable>> tri_list;
    shared_ptr<material> mat;
    Mesh(const std::vector<shared_ptr<hittable>> &_triangles) : area(0.0)
    {
        tri_list = std::vector<shared_ptr<hittable>>(_triangles);
        tri_areas.resize(_triangles.size(), 0.0);

        mesh_bvh = make_shared<bvh_node>(_triangles);
        for (size_t i = 0; i < _triangles.size(); i++)
        {
            shared_ptr<hittable> temp_tri = _triangles[i];
            area += temp_tri->getarea();
        }
        for (size_t i = 1; i < _triangles.size(); i++)
        {
            shared_ptr<hittable> temp_tri = _triangles[i];

            tri_areas[i] = tri_areas[i - 1] + temp_tri->getarea() / area;
        }
        mat = tri_list[0]->mat_ptr;
    }

    bool hasEmission()
    {
        return !(mat->getEmission().near_zero());
    }

    double getArea() const { return area; }

    bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec) const
    {
        return mesh_bvh->hit(r, t_min, t_max, rec);
    }

    bool bounding_box(double time0, double time1, aabb &output_box) const
    {
        return mesh_bvh->bounding_box(0, 0, output_box);
    }

    double pdf_value(const vec3 &o, const vec3 &v) const;
    vec3 sample(const vec3 &o) const;
};

double Mesh::pdf_value(const vec3 &o, const vec3 &v) const
{
    return 0.0;
}

vec3 Mesh::sample(const vec3 &o) const
{
    double area_ratio = random_double(0, 1.0) + EPSILON, temp_area = 0.0;
    auto it = std::lower_bound(tri_areas.begin(), tri_areas.end(), area_ratio);
    int index = std::distance(tri_areas.begin(), it) - 1;
    return tri_list[index]->sample(o);
}

#endif