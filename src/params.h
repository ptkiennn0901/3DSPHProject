/* params.h — Tham so mo phong cho kich ban Dam Break / Pool 3D.
 * Gom mot cho de chinh va de ca serial lan MPI dung chung. */
#ifndef PARAMS_H
#define PARAMS_H

#include "vec3.h"

typedef struct {
    /* Kich thuoc be (m): L = 1.0 x 0.6 x 0.6 */
    double Lx;
    double Ly;
    double Lz;

    /* Khoi nuoc ban dau */
    double waterX;
    double waterY;
    double waterZ;

    /* Roi rac hoa */
    double dx;     /* khoang cach hat */
    double h;      /* ban kinh ho tro = 2*dx */

    /* Vat ly chat luu */
    double rho0;   /* mat do nghi (kg/m^3) */
    double k;      /* do cung (he so phuong trinh trang thai) */
    double mu;     /* do nhot dong luc (Pa.s) */
    Vec3   g;      /* trong luc */

    /* Tich phan thoi gian */
    double dt;
    int    steps;
    int    out_every;

    /* Bien */
    double damping;
} Params;

void   params_default(Params *P);     /* gan gia tri mac dinh (kich ban dam) */
void   params_set_pool(Params *P);    /* nuoc lap gan day be (cho benchmark) */
void   params_finalize(Params *P);    /* dat h = 2*dx */
double params_mass(const Params *P);  /* khoi luong moi hat = rho0 * dx^3 */

#endif /* PARAMS_H */
