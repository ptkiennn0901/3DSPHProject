/* sph_core.c — Cai dat loi vat ly SPH. */
#include "sph_core.h"
#include "kernels.h"
#include <math.h>

void sph_init_block(ParticleArray *out, const Params *P)
{
    int nx = (int) (P->waterX / P->dx);
    int ny = (int) (P->waterY / P->dx);
    int nz = (int) (P->waterZ / P->dx);
    double off = 0.5 * P->dx;   /* lech nua o de hat khong nam ngay tren tuong */
    int i, j, k;

    pa_clear(out);
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                Particle pt;
                pt.pos   = vec3_make(off + i * P->dx,
                                     off + j * P->dx,
                                     off + k * P->dx);
                pt.vel   = vec3_zero();
                pt.force = vec3_zero();
                pt.rho   = 0.0;
                pt.p     = 0.0;
                pa_push(out, pt);
            }
        }
    }
}

void sph_density_pressure(ParticleArray *ps, int n_local,
                          const Grid *grid, const Params *P, IntArray *neigh)
{
    double m  = params_mass(P);
    double h2 = P->h * P->h;
    double h9 = pow(P->h, 9.0);
    int i, t;

    for (i = 0; i < n_local; i++) {
        Particle *a = &ps->data[i];
        double rho = 0.0;
        grid_neighbors(grid, a->pos, neigh);
        for (t = 0; t < neigh->size; t++) {
            int j = neigh->data[t];
            double r2 = vec3_dist2(a->pos, ps->data[j].pos);
            rho += m * kernel_poly6(r2, h2, h9);
        }
        a->rho = rho;
        a->p   = P->k * (rho - P->rho0);   /* phuong trinh trang thai */
    }
}

void sph_forces(ParticleArray *ps, int n_local,
                const Grid *grid, const Params *P, IntArray *neigh)
{
    double m  = params_mass(P);
    double h6 = pow(P->h, 6.0);
    int i, t;

    for (i = 0; i < n_local; i++) {
        Particle *a = &ps->data[i];
        Vec3 f_press = vec3_zero();
        Vec3 f_visc  = vec3_zero();
        Vec3 f_grav;

        grid_neighbors(grid, a->pos, neigh);
        for (t = 0; t < neigh->size; t++) {
            int j = neigh->data[t];
            Particle *b = &ps->data[j];
            Vec3 rij;
            double r;
            if (j == i) {
                continue;
            }
            rij = vec3_sub(a->pos, b->pos);
            r = vec3_norm(rij);
            if (r <= 0.0 || r >= P->h) {
                continue;
            }
            if (b->rho <= 1e-9) {
                continue;
            }
            /* Luc ap suat (doi xung hoa): -m (p_i+p_j)/(2 rho_j) gradW_spiky */
            {
                Vec3   grad  = kernel_spiky_grad(rij, r, P->h, h6);
                double coefp = -m * (a->p + b->p) / (2.0 * b->rho);
                f_press = vec3_add_scaled(f_press, grad, coefp);
            }
            /* Luc nhot: mu m (v_j - v_i)/rho_j lapW_visc */
            {
                double lap   = kernel_visc_lap(r, P->h, h6);
                double coefv = P->mu * m / b->rho * lap;
                Vec3   dv    = vec3_sub(b->vel, a->vel);
                f_visc = vec3_add_scaled(f_visc, dv, coefv);
            }
        }
        f_grav = vec3_scale(P->g, a->rho);   /* trong luc theo mat do */
        a->force = vec3_add(vec3_add(f_press, f_visc), f_grav);
    }
}

void sph_integrate(ParticleArray *ps, int n_local, const Params *P,
                   double wall_xmin, double wall_xmax)
{
    double d = P->damping;
    int i;

    for (i = 0; i < n_local; i++) {
        Particle *a = &ps->data[i];
        if (a->rho <= 1e-9) {
            a->rho = P->rho0;    /* an toan so hoc */
        }
        /* Mueller (2003): force la MAT DO LUC -> gia toc = F/rho. */
        a->vel = vec3_add_scaled(a->vel, a->force, P->dt / a->rho);
        a->pos = vec3_add_scaled(a->pos, a->vel, P->dt);

        /* 6 mat tuong: keo ve bien va dao dau van toc (giam damping). */
        if (a->pos.x < wall_xmin) { a->pos.x = wall_xmin; a->vel.x = -a->vel.x * d; }
        if (a->pos.x > wall_xmax) { a->pos.x = wall_xmax; a->vel.x = -a->vel.x * d; }
        if (a->pos.y < 0.0)       { a->pos.y = 0.0;       a->vel.y = -a->vel.y * d; }
        if (a->pos.y > P->Ly)     { a->pos.y = P->Ly;     a->vel.y = -a->vel.y * d; }
        if (a->pos.z < 0.0)       { a->pos.z = 0.0;       a->vel.z = -a->vel.z * d; }
        if (a->pos.z > P->Lz)     { a->pos.z = P->Lz;     a->vel.z = -a->vel.z * d; }
    }
}

double sph_density_bruteforce(const ParticleArray *ps, int i, const Params *P)
{
    double m  = params_mass(P);
    double h2 = P->h * P->h;
    double h9 = pow(P->h, 9.0);
    double rho = 0.0;
    int j;
    for (j = 0; j < ps->size; j++) {
        double r2 = vec3_dist2(ps->data[i].pos, ps->data[j].pos);
        rho += m * kernel_poly6(r2, h2, h9);
    }
    return rho;
}

double sph_self_test(ParticleArray *ps, const Grid *grid,
                     const Params *P, IntArray *neigh, int samples)
{
    double h2 = P->h * P->h;
    double h9 = pow(P->h, 9.0);
    double m  = params_mass(P);
    double max_rel = 0.0;
    int s;

    if (ps->size <= 0) {
        return 0.0;
    }
    for (s = 0; s < samples; s++) {
        int i = (s * 7919) % ps->size;   /* lay mau rai deu */
        double rho_grid = 0.0;
        double rho_bf;
        double rel;
        int t;

        grid_neighbors(grid, ps->data[i].pos, neigh);
        for (t = 0; t < neigh->size; t++) {
            int j = neigh->data[t];
            double r2 = vec3_dist2(ps->data[i].pos, ps->data[j].pos);
            rho_grid += m * kernel_poly6(r2, h2, h9);
        }
        rho_bf = sph_density_bruteforce(ps, i, P);
        rel = (rho_bf > 0.0) ? fabs(rho_grid - rho_bf) / rho_bf : 0.0;
        if (rel > max_rel) {
            max_rel = rel;
        }
    }
    return max_rel;
}
