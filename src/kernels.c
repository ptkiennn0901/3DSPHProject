/* kernels.c — Cai dat cac kernel SPH 3D. */
#include "kernels.h"

double kernel_poly6(double r2, double h2, double h9)
{
    if (r2 >= h2) {
        return 0.0;
    }
    double diff = h2 - r2;
    double coef = 315.0 / (64.0 * SPH_PI * h9);
    return coef * diff * diff * diff;
}

Vec3 kernel_spiky_grad(Vec3 rij, double r, double h, double h6)
{
    if (r <= 0.0 || r >= h) {
        return vec3_zero();
    }
    double hr   = h - r;
    double coef = -45.0 / (SPH_PI * h6) * hr * hr / r;
    return vec3_scale(rij, coef);
}

double kernel_visc_lap(double r, double h, double h6)
{
    if (r >= h) {
        return 0.0;
    }
    return 45.0 / (SPH_PI * h6) * (h - r);
}
