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
#define TINYOBJLOADER_IMPLEMENTATION
#include "rtweekend.h"

#include "camera.h"
#include "color.h"
#include "material.h"
#include "scene.h"
#include "omp.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <iostream>

struct render_params
{
    std::string scene_name;
    int image_width;
    int image_height;
    int spp;
    int spl;
    std::string getimgname(bool hasgc = false)
    {
        if (hasgc)
            return scene_name + "_" + std::to_string(spp) + "_g.png";
        else
            return scene_name + "_" + std::to_string(spp) + ".png";
    }
};

void print_render_info(render_params &rp)
{
    printf("scene: %s, width x height : [ %d, %d ], spp: %d, spl: %d\n", rp.getimgname().c_str(), rp.image_width, rp.image_height, rp.spp, rp.spl);
}

bool init(std::string xml_path, camera &cam, Scene &scene, render_params &rp);

int main(int argc, char **argv)
{

    image_texture it("../../scenes/staircase/textures/Wallpaper.jpg");
    int t_width = it.get_width(), t_height = it.get_height();
    unsigned char *t_img = new unsigned char[t_width * t_height * 3];
    for (int j = 0; j < t_height; j++)
    {
        for (int i = 0; i < t_width; i++)
        {
            int index = j * t_width + i;
            color c = it.value((i + 100 + 0.5) / t_width, (j + 200 + 0.5) / t_height, vec3());
            int scale = 255;
            t_img[index * 3] = scale * c.x();
            t_img[index * 3 + 1] = scale * c.y();
            t_img[index * 3 + 2] = scale * c.z();
        }
    }
    stbi_write_png("texture.png", t_width, t_height, 3, t_img, 0);
    // return 0;
    Scene scene;
    camera cam;
    std::string xml_path;
    render_params rp;
    std::vector<std::string> xml_paths = {
        "../../scenes/cornell-box/cornell-box.xml",
        "../../scenes/staircase/staircase.xml",
        "../../scenes/veach-mis/veach-mis.xml",
        "../../scenes/bunny-box/bunny-box.xml"};
    int test_scene = 3;
    if (argc > 1)
    {
        test_scene = std::stoi(argv[1]);
    }
    xml_path = xml_paths[test_scene];

    xml_path = xml_paths[test_scene];

    if (!init(xml_path, cam, scene, rp))
    {
        printf("read failed\n");
        return 0;
    }
    std::vector<int> ssp_test = {16, 32, 64, 128, 256, 512, 1024};
    if (argc > 2)
    {
        ssp_test = {std::stoi(argv[2])};
    }
    // std::vector<int> ssp_test = {2, 4, 8};
    for (int s = 0; s < ssp_test.size(); s++)
    {
        unsigned char *img = new unsigned char[rp.image_height * rp.image_width * 3];
        rp.spp = ssp_test[s];
        print_render_info(rp);

        for (int j = rp.image_height - 1; j >= 0; --j)
        {
            std::cerr << "\rspp: " << rp.spp << " Scanlines remaining: " << j << ' ' << std::flush;
#pragma omp parallel for
            for (int i = 0; i < rp.image_width; ++i)
            {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < rp.spp; ++s)
                {
                    // std::cout << "nmsl";
                    auto u = (i + random_double()) / (rp.image_width);
                    auto v = (j + random_double()) / (rp.image_height);
                    // auto u = float(i + 0.5) / (image_width);
                    // auto v = float(j + 0.5) / (image_height);

                    ray r = cam.get_ray(u, v);
                    // pixel_color += ray_color(r, background, scene, lights, max_depth);
                    // if (temp_c[0] < 0 || temp_c[1] < 0 || temp_c[2] < 0)
                    // {
                    //     std::cout << "nmsl";
                    // }
                    // pixel_color += scene.shade(r, 0);
                    pixel_color += scene.draw(r);
                    // if (temp_c[1] == 1.0 && temp_c[2] == 1.0)
                    // {
                    //     std::cout << "light" << pixel_color[0] << pixel_color[1] << pixel_color[2] << std::endl;
                    // }
                }
                // write_color(pixel_color, samples_per_pixel);
                write_color_xy(img, i, j, rp.image_width, pixel_color, rp.spp);
            }
        }
        stbi_flip_vertically_on_write(true);
        stbi_write_png(rp.getimgname().c_str(), rp.image_width, rp.image_height, 3, img, 0);
        // stbi_write_png(rp.getimgname(true).c_str(), rp.image_width, rp.image_height, 3, img_gamma, 0);
    }
    std::cerr << "/nDone./n";
}

bool init(std::string xml_path, camera &cam, Scene &scene, render_params &rp)
{
    tinyxml2::XMLDocument doc;
    doc.LoadFile(xml_path.c_str());
    if (doc.Error())
    {
        std::cout << "Error loading file! " << xml_path.c_str() << std::endl;
        return false;
    }
    std::string obj_path = xml_path.substr(0, xml_path.find_last_of('.')) + ".obj";

    int image_width, image_height;
    vec3 lookat, lookfrom, vup;
    float vfov, aspect, foucust_dist;
    tinyxml2::XMLElement *cam_elem = doc.FirstChildElement();
    image_width = std::stoi(cam_elem->FindAttribute("width")->Value());
    image_height = std::stoi(cam_elem->FindAttribute("height")->Value());
    vfov = std::stof(cam_elem->FindAttribute("fovy")->Value());
    for (tinyxml2::XMLElement *elem = cam_elem->FirstChildElement(); elem != NULL; elem = elem->NextSiblingElement())
    {
        const char *elemType = elem->Value();
        std::cout << "Element type: " << elemType << std::endl;
        if (strcmp(elemType, "eye") == 0)
        {
            lookfrom = vec3(std::stof(elem->FindAttribute("x")->Value()),
                            std::stof(elem->FindAttribute("y")->Value()),
                            std::stof(elem->FindAttribute("z")->Value()));
        }
        if (strcmp(elemType, "lookat") == 0)
        {
            lookat = vec3(std::stof(elem->FindAttribute("x")->Value()),
                          std::stof(elem->FindAttribute("y")->Value()),
                          std::stof(elem->FindAttribute("z")->Value()));
        }
        if (strcmp(elemType, "up") == 0)
        {
            vup = vec3(std::stof(elem->FindAttribute("x")->Value()),
                       std::stof(elem->FindAttribute("y")->Value()),
                       std::stof(elem->FindAttribute("z")->Value()));
        }
    }
    tinyxml2::XMLElement *render_elem = cam_elem->NextSiblingElement();
    while (strcmp(render_elem->Name(), "render") != 0)
    {
        render_elem = render_elem->NextSiblingElement();
    }
    std::string obj_name = xml_path.substr(xml_path.find_last_of('/') + 1);
    obj_name = obj_name.substr(0, obj_name.find_last_of('.'));
    rp.image_width = image_width;
    rp.image_height = image_height;
    rp.spp = std::stoi(render_elem->FindAttribute("spp")->Value());
    rp.spl = std::stoi(render_elem->FindAttribute("spl")->Value());
    rp.scene_name = obj_name;
    aspect = float(image_width) / image_height;
    foucust_dist = 10.f;
    cam = camera(lookfrom, lookat, vup, vfov, aspect, foucust_dist);
    scene.loadobj(obj_path);
    return true;
}
