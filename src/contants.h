#ifndef CONSTANT_H
#define CONSTANT_H

static constexpr double PI = 3.14159265358979323846;
static constexpr double INV_PI = 0.31830988618379067154;
static constexpr double HALF_PI = 1.57079632679489661923;
static constexpr double TWO_PI = 6.283185307179586476925;
static constexpr double EPSILON = 1e-8;
static constexpr int SAMPLES_PER_LIGHT = 1;
static constexpr int SAMPLES_PER_PIXEL = 64;
static constexpr double T_MAX = 1.7976931348623157e+10;
static constexpr int SHADE_DEPTH = 2;
static constexpr int BVH_MAX = 100;
static constexpr double RUSSIAN_ROUETTE = 0.95;
static constexpr double T_MIN = EPSILON;
#endif