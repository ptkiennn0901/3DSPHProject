/* vec3.h — Toan vector 3D cho mo phong SPH (ngon ngu C).
 * C khong co nap chong toan tu nen moi phep tinh la mot ham rieng. */
#ifndef VEC3_H
#define VEC3_H

typedef struct {
    double x;
    double y;
    double z;
} Vec3;

Vec3   vec3_make(double x, double y, double z);
Vec3   vec3_zero(void);
Vec3   vec3_add(Vec3 a, Vec3 b);
Vec3   vec3_sub(Vec3 a, Vec3 b);
Vec3   vec3_scale(Vec3 a, double s);
Vec3   vec3_div(Vec3 a, double s);
Vec3   vec3_negate(Vec3 a);
Vec3   vec3_add_scaled(Vec3 a, Vec3 b, double s); /* a + b*s */
double vec3_dot(Vec3 a, Vec3 b);
Vec3   vec3_cross(Vec3 a, Vec3 b);
double vec3_norm2(Vec3 a);
double vec3_norm(Vec3 a);
double vec3_dist2(Vec3 a, Vec3 b);
double vec3_dist(Vec3 a, Vec3 b);
Vec3   vec3_normalize(Vec3 a);
Vec3   vec3_lerp(Vec3 a, Vec3 b, double t);     /* noi suy tuyen tinh */
Vec3   vec3_min(Vec3 a, Vec3 b);
Vec3   vec3_max(Vec3 a, Vec3 b);
Vec3   vec3_abs(Vec3 a);
int    vec3_inside_box(Vec3 p, Vec3 lo, Vec3 hi);

#endif /* VEC3_H */
