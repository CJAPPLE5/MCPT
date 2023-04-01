#ifndef Scene_H
#define Scene_H
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

#include "rtweekend.h"

#include "hittable.h"
#include "triangle.h"
#include "bvh.h"
#include "texture.h"
#include "material.h"
#include "mesh.h"

#include <memory>
#include <vector>

double move(double tex_raw)
{
    if (tex_raw < 0)
    {
        tex_raw = -tex_raw;
        tex_raw = std::floor(tex_raw) + std::ceil(tex_raw) - tex_raw;
    }
    if (tex_raw > 1)
    {
        tex_raw = tex_raw - std::floor(tex_raw);
    }
    return tex_raw;
}

class Scene : public hittable
{
public:
    Scene()
    {
        mesh_count = 0;
        environment = make_shared<solid_color>(0.2, 0.2, 0.2);
    };

    bool loadobj(std::string obj_path);

    void clear()
    {
        materials.clear();
        meshes.clear();
    }

    color shade(const ray &r, int depth);
    color shade_(const ray &r, hit_record &rec, int depth);
    color draw(const ray &r);
    bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec) const;

    bool bounding_box(double time0, double time1, aabb &output_box) const;
    double pdf_value(const vec3 &o, const vec3 &v) const;
    vec3 sample(const vec3 &o) const;
    void sample_light(vec3 &light_position, double &pdf, double &light_area);

private:
    // std::vector<shared_ptr<hittable>> elements;
    // std::vector<shared_ptr<hittable>> elements;
    std::map<std::string, shared_ptr<material>> materials;
    std::vector<shared_ptr<Mesh>> meshes;
    std::vector<shared_ptr<Mesh>> lights;
    std::vector<double> lights_area_ratio;
    // shared_ptr<bvh_node> bvh_root;
    shared_ptr<texture> environment;
    double lights_area;
    int mesh_count;
};

void Scene::sample_light(vec3 &light_position, double &pdf, double &light_area)
{
    double p = random_double();
    auto iter = std::lower_bound(lights_area_ratio.begin(), lights_area_ratio.end(), p);
    if (iter == lights_area_ratio.end())
    {
        iter = lights_area_ratio.begin();
    }
    int index = std::distance(lights_area_ratio.begin(), iter);
    light_position = lights[index]->sample(vec3());
    pdf = lights[index]->getArea() / lights_area;
    light_area = lights[index]->getArea();
}

color Scene::draw(const ray &r)
{
    double t_min = T_MIN, t_max = T_MAX;
    hit_record hit_rec;
    if (hit(r, t_min, t_max, hit_rec))
    {
        return shade_(r, hit_rec, 0);
    }
    else
    {
        return color();
    }
}

color Scene::shade(const ray &r, int depth)
{
    color l_total, l_e, l_direct, l_indirect;
    hit_record rec;
    if (depth == MAX_DEPTH)
    {
        return color();
    }
    if (!hit(r, T_MIN, T_MAX, rec))
    {
        return color();
    }

    if (depth == 0 && rec.mat_ptr->hasemission())
    {
        l_e = rec.mat_ptr->getEmission();
    }
    if (rec.mat_ptr->hasTexture())
    {
        return rec.mat_ptr->map_kd->value(rec.uv.x(), rec.uv.y(), rec.p);
    }

    if (!rec.mat_ptr->hasRefraction())
    {
        rec.normal = dot(rec.normal, -r.dir) >= 0 ? rec.normal : -rec.normal;
    }

    // {
    //     double light_pdf, light_area;
    //     point3 light_position;
    //     sample_light(light_position, light_pdf, light_area);

    //     vec3 light_dir = normalize(light_position - rec.p);
    //     hit_record temp_rec;
    // ray shadow_ray(rec.p, light_dir);
    //     if (hit(shadow_ray, T_MIN, T_MAX, temp_rec) && temp_rec.mat_ptr->hasemission() &&
    //         dot(temp_rec.normal, light_dir) < 0 && (temp_rec.p - light_position).length() < EPSILON)
    //     {
    //         double cos0 = fmax(dot(rec.normal, light_dir), 0.0);
    //         double cos1 = fmax(dot(temp_rec.normal, -light_dir), 0.0);
    //         vec3 brdf;
    //         scatter_record temp_srec;
    //         temp_srec.scatter_ray.orig = rec.p;
    //         temp_srec.scatter_ray.dir = light_dir;
    //         brdf = rec.mat_ptr->brdf(r, rec, temp_srec);
    //         vec3 direct;
    //         if (cos0 > 0 && cos1 > 0)
    //         {
    //             direct = brdf * temp_rec.mat_ptr->getEmission() * cos0 * cos1 *
    //                      light_area / (light_position - rec.p).length_squared();
    //             l_direct += direct;
    //         }
    //     }
    // }
    if (depth != 0)
        for (auto &mesh : meshes)
        {
            color direct_light;
            if (mesh->hasEmission())
            {
                for (size_t i = 0; i < SAMPLES_PER_LIGHT; i++)
                {
                    point3 light_position = mesh->sample(rec.normal);

                    vec3 light_dir = normalize(light_position - rec.p);
                    hit_record temp_rec;
                    ray shadow_ray(rec.p, light_dir);
                    if (hit(shadow_ray, T_MIN, T_MAX, temp_rec) && temp_rec.mesh_id == mesh->mesh_id &&
                        dot(temp_rec.normal, light_dir) < 0 && (temp_rec.p - light_position).length() < EPSILON)
                    {
                        double cos0 = fmax(dot(rec.normal, light_dir), 0.0);
                        double cos1 = fmax(dot(temp_rec.normal, -light_dir), 0.0);
                        vec3 brdf;
                        scatter_record temp_srec;
                        temp_srec.scatter_ray.orig = rec.p;
                        temp_srec.scatter_ray.dir = light_dir;
                        brdf = rec.mat_ptr->brdf(r, rec, temp_srec);
                        vec3 direct;
                        if (cos0 > 0 && cos1 > 0)
                        {
                            direct = brdf * temp_rec.mat_ptr->getEmission() * cos0 * cos1 * mesh->getArea() / (light_position - rec.p).length_squared();
                            direct_light += direct;
                        }
                    }
                }
                l_direct += direct_light / SAMPLES_PER_LIGHT;
            }
        }

    // return l_direct + l_e;
    if (random_double() < RUSSIAN_ROUETTE)
    // if (0)
    {
        scatter_record srec;
        rec.mat_ptr->scatter_(r, rec, srec);
        double cos0 = dot(rec.normal, srec.scatter_ray.dir);
        // hit_record rec_temp;
        // hit(srec.scatter_ray, T_MIN, T_MAX, rec_temp);

        if (srec.scatter_type == REFRACT)
        {
            l_indirect += shade(srec.scatter_ray, depth) * rec.mat_ptr->tr / srec.pdf / RUSSIAN_ROUETTE;
        }
        else
        {
            if (cos0 > 0 && srec.pdf > 0)
            {
                vec3 coeff = rec.mat_ptr->brdf(r, rec, srec) * cos0 / srec.pdf / RUSSIAN_ROUETTE;
                l_indirect += coeff * shade(srec.scatter_ray, depth + 1);
                // l_indirect += coeff * rec.mat_ptr->kd;
            }
        }
    }
    l_total = l_e + l_direct + l_indirect;
    return l_total;
}

color Scene::shade_(const ray &r, hit_record &rec, int depth)
{
    color l_total, l_e, l_direct, l_indirect;

    if (rec.mat_ptr->hasemission())
    {
        return rec.mat_ptr->getEmission();
    }
    // if (!rec.mat_ptr->hasRefraction())
    // {
    //     rec.normal = dot(rec.normal, -r.dir) >= 0 ? rec.normal : -rec.normal;
    // }

    for (auto &mesh : meshes)
    {
        color direct_light;
        if (mesh->hasEmission())
        {
            for (size_t i = 0; i < SAMPLES_PER_LIGHT; i++)
            {
                point3 light_position = mesh->sample(rec.normal);

                vec3 light_dir = normalize(light_position - rec.p);
                hit_record temp_rec;
                ray shadow_ray(rec.p, light_dir);
                if (hit(shadow_ray, T_MIN, T_MAX, temp_rec) && temp_rec.mesh_id == mesh->mesh_id &&
                    dot(temp_rec.normal, light_dir) < 0 && (temp_rec.p - light_position).length() < EPSILON)
                {
                    double cos0 = fmax(dot(rec.normal, light_dir), 0.0);
                    double cos1 = fmax(dot(temp_rec.normal, -light_dir), 0.0);
                    vec3 brdf;
                    scatter_record temp_srec;
                    temp_srec.scatter_ray.orig = rec.p;
                    temp_srec.scatter_ray.dir = light_dir;
                    brdf = rec.mat_ptr->brdf(r, rec, temp_srec);
                    vec3 direct;
                    if (cos0 > 0 && cos1 > 0)
                    {
                        direct = brdf * temp_rec.mat_ptr->getEmission() * cos0 * cos1 * mesh->getArea() / (light_position - rec.p).length_squared();
                        direct_light += direct;
                    }
                }
            }
            l_direct += direct_light / SAMPLES_PER_LIGHT;
        }
    }

    // return l_direct + l_e;
    if (random_double() < RUSSIAN_ROUETTE)
    // if (0)
    {
        scatter_record srec;
        rec.mat_ptr->scatter_(r, rec, srec);
        double cos0 = dot(rec.normal, srec.scatter_ray.dir);
        // hit_record rec_temp;
        // hit(srec.scatter_ray, T_MIN, T_MAX, rec_temp);
        if (srec.pdf > 0)
        {
            hit_record temp_rec;
            if (hit(srec.scatter_ray, T_MIN, T_MAX, temp_rec) && !temp_rec.mat_ptr->hasemission())
            {
                if (srec.scatter_type == REFRACT)
                {
                    l_indirect += shade(srec.scatter_ray, depth) * rec.mat_ptr->tr / srec.pdf / RUSSIAN_ROUETTE;
                }
                else
                {
                    if (cos0 > 0 && srec.pdf > 0)
                    {
                        vec3 coeff = rec.mat_ptr->brdf(r, rec, srec) * cos0 / srec.pdf / RUSSIAN_ROUETTE;
                        l_indirect += coeff * shade(srec.scatter_ray, depth + 1);
                        // l_indirect += coeff * rec.mat_ptr->kd;
                    }
                }
            }
        }
    }
    l_total = l_e + l_direct + l_indirect;
    return l_total;
}

bool Scene::loadobj(std::string obj_path)
{
    lights_area = 0.0;
    std::cout << "read obj:" << obj_path << std::endl;
    std::map<std::string, vec3> light_radinaces;
    std::string base_path = obj_path.substr(0, obj_path.find_last_of('/'));
    std::string scene_name = obj_path.substr(obj_path.find_last_of('/') + 1);
    std::string xml_path = base_path + "/" + scene_name.substr(0, scene_name.find_last_of('.')) + ".xml";
    tinyxml2::XMLDocument xml_doc;
    xml_doc.LoadFile(xml_path.c_str());
    if (xml_doc.Error())
    {
        std::cout << "Error loading file!" << xml_path.c_str() << std::endl;
        return 1;
    }

    for (tinyxml2::XMLElement *lightElem = xml_doc.FirstChildElement("light"); lightElem; lightElem = lightElem->NextSiblingElement("light"))
    {
        const char *mtlname = lightElem->Attribute("mtlname");
        const char *radiance = lightElem->Attribute("radiance");

        std::vector<double> radiance_values;

        std::stringstream ss(radiance);
        std::string item;
        while (std::getline(ss, item, ','))
        {
            radiance_values.push_back(std::stof(item));
        }
        light_radinaces[std::string(mtlname)] = vec3(radiance_values[0], radiance_values[1], radiance_values[2]);
    }

    // tinyxml2::XMLElement *root = xml_doc.RootElement();
    // for (tinyxml2::XMLElement *element = root->FirstChildElement("light"); element != NULL;
    //      element = element->NextSiblingElement())
    // {
    //     const char *light_name = element->Attribute("mtlname");
    //     const char *light_radinace_c = element->Attribute("radiance");
    //     std::vector<std::string> values;
    //     std::vector<double> radinace_values;
    //     std::stringstream ss(light_radinace_c);
    //     std::string item;
    //     while (std::getline(ss, item, ','))
    //     {
    //         radinace_values.push_back(std::stoi(item));
    //     }
    //     light_radinaces[std::string(light_name)] = vec3(radinace_values[0], radinace_values[1], radinace_values[2]);
    // }
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> obj_shapes;
    std::vector<tinyobj::material_t> obj_materials;

    std::string warn;
    std::string err;
    bool ret = tinyobj::LoadObj(&attrib, &obj_shapes, &obj_materials, &warn, &err, obj_path.c_str(),
                                base_path.c_str(), true);

    // if (!warn.empty())
    // {
    //     std::cout << "WARN: " << warn << std::endl;
    // }
    // if (!err.empty())
    // {
    //     std::cerr << "ERR: " << err << std::endl;
    // }
    if (!ret)
    {
        printf("Failed to load/parse .obj.\n");
        return false;
    }

    for (size_t i = 0; i < obj_materials.size(); i++)
    {
        vec3 kd(static_cast<const double>(obj_materials[i].diffuse[0]),
                static_cast<const double>(obj_materials[i].diffuse[1]),
                static_cast<const double>(obj_materials[i].diffuse[2])),
            ks(static_cast<const double>(obj_materials[i].specular[0]),
               static_cast<const double>(obj_materials[i].specular[1]),
               static_cast<const double>(obj_materials[i].specular[2])),
            tr(static_cast<const double>(obj_materials[i].transmittance[0]),
               static_cast<const double>(obj_materials[i].transmittance[1]),
               static_cast<const double>(obj_materials[i].transmittance[2]));
        double ns = static_cast<const double>(obj_materials[i].shininess),
               ni = static_cast<const double>(obj_materials[i].ior);
        vec3 ke;
        if (light_radinaces.find(obj_materials[i].name) != light_radinaces.end())
        {
            ke = light_radinaces[obj_materials[i].name];
        }

        shared_ptr<image_texture> temp_texture;
        bool has_texture = false;
        if (obj_materials[i].diffuse_texname.length() > 0)
        {
            std::cout << "read texture:" << obj_materials[i].diffuse_texname << std::endl;
            temp_texture = make_shared<image_texture>((base_path + "/" + obj_materials[i].diffuse_texname).c_str());
            has_texture = true;
        }
        else
        {
            temp_texture = nullptr;
        }
        materials[obj_materials[i].name] = make_shared<material>(kd, ks, ke, ni, ns, tr, has_texture, temp_texture);
    }

    // meshes.resize(obj_shapes.size());

    for (size_t i = 0; i < obj_shapes.size(); i++)
    {
        size_t index_offset = 0;
        std::vector<shared_ptr<hittable>> temp_triangles;
        temp_triangles.resize(obj_shapes[i].mesh.num_face_vertices.size());
        // size_t base_size = elements.size();
        // elements.resize(base_size + obj_shapes[i].mesh.num_face_vertices.size());
        // meshes[i].resize(obj_shapes[i].mesh.num_face_vertices.size());
        // face
        // if (obj_shapes[i].mesh.num_face_vertices.size() > 100)
        // {
        //     continue;
        // }
        for (size_t f = 0; f < obj_shapes[i].mesh.num_face_vertices.size(); f++)
        {
            size_t fnum = obj_shapes[i].mesh.num_face_vertices[f];
            // assert(fnum == 3, "fnum != 3");
            tinyobj::index_t idx0 = obj_shapes[i].mesh.indices[index_offset + 0];
            tinyobj::index_t idx1 = obj_shapes[i].mesh.indices[index_offset + 1];
            tinyobj::index_t idx2 = obj_shapes[i].mesh.indices[index_offset + 2];
            point3 p0(static_cast<double>(attrib.vertices[idx0.vertex_index * 3]),
                      static_cast<double>(attrib.vertices[idx0.vertex_index * 3 + 1]),
                      static_cast<double>(attrib.vertices[idx0.vertex_index * 3 + 2])),
                p1(static_cast<double>(attrib.vertices[idx1.vertex_index * 3]),
                   static_cast<double>(attrib.vertices[idx1.vertex_index * 3 + 1]),
                   static_cast<double>(attrib.vertices[idx1.vertex_index * 3 + 2])),
                p2(static_cast<double>(attrib.vertices[idx2.vertex_index * 3]),
                   static_cast<double>(attrib.vertices[idx2.vertex_index * 3 + 1]),
                   static_cast<double>(attrib.vertices[idx2.vertex_index * 3 + 2]));

            normal3 n0(static_cast<double>(attrib.normals[idx0.normal_index * 3]),
                       static_cast<double>(attrib.normals[idx0.normal_index * 3 + 1]),
                       static_cast<double>(attrib.normals[idx0.normal_index * 3 + 2])),
                n1(static_cast<double>(attrib.normals[idx1.normal_index * 3]),
                   static_cast<double>(attrib.normals[idx1.normal_index * 3 + 1]),
                   static_cast<double>(attrib.normals[idx1.normal_index * 3 + 2])),
                n2(static_cast<double>(attrib.normals[idx2.normal_index * 3]),
                   static_cast<double>(attrib.normals[idx2.normal_index * 3 + 1]),
                   static_cast<double>(attrib.normals[idx2.normal_index * 3 + 2]));

            vec3 tex0(static_cast<double>(attrib.texcoords[idx0.texcoord_index * 2]),
                      static_cast<double>(attrib.texcoords[idx0.texcoord_index * 2 + 1]), 0.0);
            vec3 tex1(static_cast<double>(attrib.texcoords[idx1.texcoord_index * 2]),
                      static_cast<double>(attrib.texcoords[idx1.texcoord_index * 2 + 1]), 0.0);
            vec3 tex2(static_cast<double>(attrib.texcoords[idx2.texcoord_index * 2]),
                      static_cast<double>(attrib.texcoords[idx2.texcoord_index * 2 + 1]), 0.0);

            std::string material_name = obj_materials[obj_shapes[i].mesh.material_ids[f]].name;
            // elements[base_size + index_offset / 3] = make_shared<triangle>(p0, p1, p2, n0, n1, n2, materials[material_name]);
            // if (materials.find(material_name) == materials.end())
            // {
            //     int a = 2;
            //     temp_triangles[f] = make_shared<hittable>(p0, p1, p2, n0, n1, n2, materials[material_name]);
            // }
            temp_triangles[f] = make_shared<hittable>(p0, p1, p2, n0, n1, n2, tex0, tex1, tex2, materials[material_name]);
            index_offset += fnum;
        }
        shared_ptr<Mesh> temp_mesh = make_shared<Mesh>(temp_triangles, mesh_count++);
        meshes.push_back(temp_mesh);
        if (temp_mesh->hasEmission())
        {
            lights.push_back(temp_mesh);
            if (lights_area_ratio.empty())
            {
                lights_area_ratio.push_back(temp_mesh->getArea());
            }
            else
            {
                lights_area_ratio.push_back(lights_area_ratio.back() + temp_mesh->getArea());
            }
            lights_area += temp_mesh->getArea();
        }
    }
    for (int i = 1; i < lights_area_ratio.size(); i++)
    {
        lights_area_ratio[i] /= lights_area;
    }
    return true;
}

bool Scene::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;
    for (int i = 0; i < meshes.size(); i++)
    {
        if (meshes[i]->hit(r, t_min, closest_so_far, temp_rec))
        {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

bool Scene::bounding_box(double time0, double time1, aabb &output_box) const
{
    if (meshes.empty())
        return false;

    aabb temp_box;
    bool first_box = true;

    for (const auto &mesh : meshes)
    {
        if (!mesh->bounding_box(time0, time1, temp_box))
            return false;
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }

    return true;
}

double Scene::pdf_value(const point3 &o, const vec3 &v) const
{
    auto weight = 1.0 / meshes.size();
    auto sum = 0.0;

    for (const auto &mesh : meshes)
        sum += weight * mesh->pdf_value(o, v);

    return sum;
}

vec3 Scene::sample(const vec3 &o) const
{
    auto int_size = static_cast<int>(meshes.size());
    return meshes[random_int(0, int_size - 1)]->sample(o);
}

#endif
