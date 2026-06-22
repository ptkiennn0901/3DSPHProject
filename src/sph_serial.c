/* sph_serial.c — Ban SERIAL (1 tien trinh) cua SPH 3D.
 * La nen tang, phai dung truoc khi lam MPI.
 *
 * Dung: ./sph_serial [steps] [dx] [dt] [scenario]
 *   steps    : so buoc thoi gian
 *   dx       : khoang cach hat (nho hon -> nhieu hat hon)
 *   dt       : buoc thoi gian (giam khi dx nho: dt ~ 0.01*dx)
 *   scenario : "dam" (mac dinh, khoi goc) hoac "pool" (day be) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "params.h"
#include "darray.h"
#include "grid.h"
#include "sph_core.h"
#include "io.h"
#include "stats.h"

int main(int argc, char **argv)
{
    Params P;
    ParticleArray ps;
    IntArray neigh;
    Grid grid;
    int step;
    int frame = 0;
    int use_pool = 0;
    struct timespec t0, t1;
    double sec;

    params_default(&P);
    if (argc >= 2) {
        P.steps = atoi(argv[1]);
    }
    if (argc >= 3) {
        P.dx = atof(argv[2]);
    }
    if (argc >= 4) {
        P.dt = atof(argv[3]);
    }
    if (argc >= 5 && strcmp(argv[4], "pool") == 0) {
        params_set_pool(&P);
        use_pool = 1;
    }
    params_finalize(&P);

    pa_init(&ps);
    ia_init(&neigh);
    sph_init_block(&ps, &P);

    printf("[serial] N = %d hat | dx=%.3f h=%.3f | %d buoc | dt=%.5f | scenario=%s\n",
           ps.size, P.dx, P.h, P.steps, P.dt, use_pool ? "pool" : "dam");

    grid_init(&grid, 1 << 16, P.h, vec3_zero());

    /* Tu kiem tra: mat do tinh bang luoi phai trung mat do O(N^2). */
    grid_build(&grid, ps.data, ps.size);
    {
        double err = sph_self_test(&ps, &grid, &P, &neigh, 20);
        printf("[serial] Self-test luoi vs O(N^2): sai so lon nhat = %.2e %s\n",
               err, err < 1e-9 ? "(PASS)" : "(KIEM TRA LAI)");
    }

    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (step = 0; step < P.steps; step++) {
        grid_build(&grid, ps.data, ps.size);
        sph_density_pressure(&ps, ps.size, &grid, &P, &neigh);
        sph_forces(&ps, ps.size, &grid, &P, &neigh);
        sph_integrate(&ps, ps.size, &P, 0.0, P.Lx);

        if (step % P.out_every == 0) {
            Stats st;
            char tag[32];
            io_write_frame("output", frame, ps.data, ps.size, 0.0, 1);
            io_write_vtk("output", frame, ps.data, ps.size);
            stats_compute(&st, ps.data, ps.size, params_mass(&P));
            snprintf(tag, sizeof(tag), "serial f%03d", frame);
            stats_print(&st, tag);
            frame++;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    sec = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    printf("[serial] Xong %d buoc trong %.3f s (%.1f buoc/s)\n",
           P.steps, sec, P.steps / sec);

    grid_free(&grid);
    ia_free(&neigh);
    pa_free(&ps);
    return 0;
}
