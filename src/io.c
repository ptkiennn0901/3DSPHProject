/* io.c — Cai dat ghi file CSV. */
#include "io.h"
#include <stdio.h>
#include <math.h>

void io_write_frame(const char *dir, int frame,
                    const Particle *ps, int n, double own_axis_width, int nproc)
{
    char path[512];
    FILE *f;
    int i;

    snprintf(path, sizeof(path), "%s/frame_%04d.csv", dir, frame);
    f = fopen(path, "w");
    if (f == NULL) {
        return;
    }
    fprintf(f, "x,y,z,vx,vy,vz,rho,p,rank\n");
    for (i = 0; i < n; i++) {
        const Particle *a = &ps[i];
        int rank = 0;
        if (own_axis_width > 0.0) {
            rank = (int) floor(a->pos.x / own_axis_width);
            if (rank < 0) {
                rank = 0;
            }
            if (rank >= nproc) {
                rank = nproc - 1;
            }
        }
        fprintf(f, "%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.2f,%.2f,%d\n",
                a->pos.x, a->pos.y, a->pos.z,
                a->vel.x, a->vel.y, a->vel.z, a->rho, a->p, rank);
    }
    fclose(f);
}

void io_write_timing(const char *path, const double *compute,
                     const double *comm, int nproc)
{
    FILE *f = fopen(path, "w");
    int r;
    if (f == NULL) {
        return;
    }
    fprintf(f, "rank,compute_s,comm_s\n");
    for (r = 0; r < nproc; r++) {
        fprintf(f, "%d,%.4f,%.4f\n", r, compute[r], comm[r]);
    }
    fclose(f);
}

void io_write_vtk(const char *dir, int frame, const Particle *ps, int n)
{
    char path[512];
    FILE *f;
    int i;

    snprintf(path, sizeof(path), "%s/frame_%04d.vtk", dir, frame);
    f = fopen(path, "w");
    if (f == NULL) {
        return;
    }
    /* Header VTK legacy. */
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "SPH 3D frame %d\n", frame);
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET POLYDATA\n");

    /* Toa do diem. */
    fprintf(f, "POINTS %d float\n", n);
    for (i = 0; i < n; i++) {
        fprintf(f, "%.5f %.5f %.5f\n", ps[i].pos.x, ps[i].pos.y, ps[i].pos.z);
    }

    /* Moi diem la mot vertex. */
    fprintf(f, "VERTICES %d %d\n", n, 2 * n);
    for (i = 0; i < n; i++) {
        fprintf(f, "1 %d\n", i);
    }

    /* Truong vo huong: toc do va mat do. */
    fprintf(f, "POINT_DATA %d\n", n);
    fprintf(f, "SCALARS speed float 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (i = 0; i < n; i++) {
        double sp = vec3_norm(ps[i].vel);
        fprintf(f, "%.5f\n", sp);
    }
    fprintf(f, "SCALARS density float 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (i = 0; i < n; i++) {
        fprintf(f, "%.3f\n", ps[i].rho);
    }
    fclose(f);
}
