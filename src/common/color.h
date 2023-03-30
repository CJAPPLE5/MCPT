#ifndef COLOR_H
#define COLOR_H
//==============================================================================================
// Originally written in 2020 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "vec3.h"

#include <iostream>
#include "math.h"
void write_color_xy(unsigned char *img, int x, int y, int width, color pixel_color, int samples_pre_pixel)
{
    int index = (y * width + x) * 3;
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Replace NaN components with zero. See explanation in Ray Tracing: The Rest of Your Life.
    if (r != r)
        r = 0.0;
    if (g != g)
        g = 0.0;
    if (b != b)
        b = 0.0;

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_pre_pixel;
    r = scale * r;
    g = scale * g;
    b = scale * b;
    r = std::powf(r, 1 / 2.2);
    g = std::powf(g, 1 / 2.2);
    b = std::powf(b, 1 / 2.2);
    img[index] = 255 * clamp(r, 0.0, 0.999999999);
    img[index + 1] = 255 * clamp(g, 0.0, 0.999999999);
    img[index + 2] = 255 * clamp(b, 0.0, 0.999999999);
}

void write_color(color pixel_color, int samples_per_pixel)
{
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Replace NaN components with zero. See explanation in Ray Tracing: The Rest of Your Life.
    if (r != r)
        r = 0.0;
    if (g != g)
        g = 0.0;
    if (b != b)
        b = 0.0;

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated [0,255] value of each color component.
    std::cout << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
              << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
              << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

#endif
