/* params.c — Gia tri mac dinh va ham tien ich cho tham so mo phong. */
#include "params.h"

void params_default(Params *P)
{
    P->Lx = 1.0;
    P->Ly = 0.6;
    P->Lz = 0.6;

    /* Kich ban "dam": khoi nuoc o goc trai be (minh hoa vat ly). */
    P->waterX = 0.3;
    P->waterY = 0.4;
    P->waterZ = 0.4;

    P->dx = 0.02;
    P->h  = 0.04;

    P->rho0 = 1000.0;
    P->k    = 2000.0;
    P->mu   = 1.0;
    P->g    = vec3_make(0.0, 0.0, -9.81);

    P->dt        = 0.0004;
    P->steps     = 2000;
    P->out_every = 20;

    P->damping = 0.5;
}

void params_set_pool(Params *P)
{
    /* Nuoc lap gan day be -> hat phan bo deu theo truc x -> can bang tai. */
    P->waterX = 0.95;
    P->waterY = 0.55;
    P->waterZ = 0.50;
}

void params_finalize(Params *P)
{
    P->h = 2.0 * P->dx;
}

double params_mass(const Params *P)
{
    return P->rho0 * P->dx * P->dx * P->dx;
}
