/* particle.h — Cau truc mot hat SPH.
 * La POD (chi gom double lien tiep, khong con tro) nen gui qua MPI bang MPI_BYTE an toan. */
#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec3.h"

typedef struct {
    Vec3   pos;    /* vi tri (x,y,z) */
    Vec3   vel;    /* van toc */
    Vec3   force;  /* luc tich luy trong moi buoc */
    double rho;    /* mat do */
    double p;      /* ap suat */
} Particle;

#endif /* PARTICLE_H */
