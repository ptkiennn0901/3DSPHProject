/* kernels.h — Kernel SPH 3D theo Mueller, Charypar, Gross (2003).
 * Hang so chuan hoa cho khong gian 3 chieu.
 *   poly6      -> mat do
 *   spiky grad -> luc ap suat
 *   visc lap   -> luc nhot */
#ifndef KERNELS_H
#define KERNELS_H

#include "vec3.h"

#define SPH_PI 3.14159265358979323846

/* W_poly6(r,h) = 315/(64 pi h^9) (h^2 - r^2)^3, voi 0 <= r <= h.
 * Nhan r2 = r^2, h2 = h^2, h9 = h^9 de tranh tinh lai. */
double kernel_poly6(double r2, double h2, double h9);

/* Gradient kernel Spiky: -45/(pi h^6) (h-r)^2 * (rij/r). */
Vec3 kernel_spiky_grad(Vec3 rij, double r, double h, double h6);

/* Laplacian kernel Viscosity: 45/(pi h^6) (h-r). */
double kernel_visc_lap(double r, double h, double h6);

#endif /* KERNELS_H */
