/* vec3.c — Cai dat cac phep toan vector 3D. */
#include "vec3.h"
#include <math.h>

Vec3 vec3_make(double x, double y, double z)
{
    Vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

Vec3 vec3_zero(void)
{
    return vec3_make(0.0, 0.0, 0.0);
}

Vec3 vec3_add(Vec3 a, Vec3 b)
{
    return vec3_make(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3 vec3_sub(Vec3 a, Vec3 b)
{
    return vec3_make(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec3 vec3_scale(Vec3 a, double s)
{
    return vec3_make(a.x * s, a.y * s, a.z * s);
}

Vec3 vec3_div(Vec3 a, double s)
{
    return vec3_make(a.x / s, a.y / s, a.z / s);
}

Vec3 vec3_negate(Vec3 a)
{
    return vec3_make(-a.x, -a.y, -a.z);
}

Vec3 vec3_add_scaled(Vec3 a, Vec3 b, double s)
{
    return vec3_make(a.x + b.x * s, a.y + b.y * s, a.z + b.z * s);
}

double vec3_dot(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 vec3_cross(Vec3 a, Vec3 b)
{
    return vec3_make(a.y * b.z - a.z * b.y,
                     a.z * b.x - a.x * b.z,
                     a.x * b.y - a.y * b.x);
}

double vec3_norm2(Vec3 a)
{
    return a.x * a.x + a.y * a.y + a.z * a.z;
}

double vec3_norm(Vec3 a)
{
    return sqrt(vec3_norm2(a));
}

double vec3_dist2(Vec3 a, Vec3 b)
{
    return vec3_norm2(vec3_sub(a, b));
}

double vec3_dist(Vec3 a, Vec3 b)
{
    return sqrt(vec3_dist2(a, b));
}

Vec3 vec3_normalize(Vec3 a)
{
    double n = vec3_norm(a);
    if (n <= 1e-12) {
        return vec3_zero();
    }
    return vec3_div(a, n);
}

Vec3 vec3_lerp(Vec3 a, Vec3 b, double t)
{
    return vec3_make(a.x + (b.x - a.x) * t,
                     a.y + (b.y - a.y) * t,
                     a.z + (b.z - a.z) * t);
}

Vec3 vec3_min(Vec3 a, Vec3 b)
{
    return vec3_make(a.x < b.x ? a.x : b.x,
                     a.y < b.y ? a.y : b.y,
                     a.z < b.z ? a.z : b.z);
}

Vec3 vec3_max(Vec3 a, Vec3 b)
{
    return vec3_make(a.x > b.x ? a.x : b.x,
                     a.y > b.y ? a.y : b.y,
                     a.z > b.z ? a.z : b.z);
}

Vec3 vec3_abs(Vec3 a)
{
    return vec3_make(a.x < 0.0 ? -a.x : a.x,
                     a.y < 0.0 ? -a.y : a.y,
                     a.z < 0.0 ? -a.z : a.z);
}

int vec3_inside_box(Vec3 p, Vec3 lo, Vec3 hi)
{
    if (p.x < lo.x || p.x > hi.x) { return 0; }
    if (p.y < lo.y || p.y > hi.y) { return 0; }
    if (p.z < lo.z || p.z > hi.z) { return 0; }
    return 1;
}
