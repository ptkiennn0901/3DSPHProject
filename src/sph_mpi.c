/* sph_mpi.c — Ban MPI MULTI-NODE cua SPH 3D.
 *
 * Phan ra mien 1D theo truc x ("slab"): domain x in [0,Lx] chia P lat.
 * Rank r so huu lat [r*W, (r+1)*W) voi W = Lx/P.
 *
 * Moi buoc thoi gian:
 *   1) Ghost exchange lan 1 (vi tri)   -> tinh mat do/ap suat hat cuc bo
 *   2) Ghost exchange lan 2 (rho,p)    -> tinh luc dung cho hat bien
 *   3) Luc + tich phan
 *   4) Migration (hat vuot bien lat -> rank lan can)
 *   5) Output (MPI_Gatherv ve rank 0)
 *
 * Dung: mpirun -np P ./sph_mpi [steps] [dx] [dt] [scenario] */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "params.h"
#include "darray.h"
#include "grid.h"
#include "sph_core.h"
#include "io.h"
#include "stats.h"

/* Gui 'send' (send_n hat) toi dst, nhan tu src, NOI vao recv.
 * dst/src co the la MPI_PROC_NULL (rank bien) -> coi nhu khong gui/nhan.
 * Pha 1: trao doi count; pha 2: data (MPI_BYTE). */
static void exchange(const Particle *send, int send_n, int dst, int src,
                     ParticleArray *recv)
{
    int recv_n = 0;
    MPI_Sendrecv(&send_n, 1, MPI_INT, dst, 0,
                 &recv_n, 1, MPI_INT, src, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (recv_n > 0) {
        int old = recv->size;
        pa_resize(recv, old + recv_n);
        MPI_Sendrecv((void *) send, send_n * (int) sizeof(Particle), MPI_BYTE, dst, 1,
                     recv->data + old, recv_n * (int) sizeof(Particle), MPI_BYTE, src, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        /* Van phai goi de bat tay neu phia kia gui. */
        MPI_Sendrecv((void *) send, send_n * (int) sizeof(Particle), MPI_BYTE, dst, 1,
                     NULL, 0, MPI_BYTE, src, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/* Xay danh sach halo: hat cuc bo sat bien trai/phai. */
static void build_halo(const ParticleArray *ps, int n_local,
                       double s_min, double s_max, double h,
                       ParticleArray *to_left, ParticleArray *to_right)
{
    int i;
    pa_clear(to_left);
    pa_clear(to_right);
    for (i = 0; i < n_local; i++) {
        double x = ps->data[i].pos.x;
        if (x < s_min + h) {
            pa_push(to_left, ps->data[i]);
        }
        if (x > s_max - h) {
            pa_push(to_right, ps->data[i]);
        }
    }
}

int main(int argc, char **argv)
{
    int rank;
    int P;
    Params par;
    int use_pool = 0;

    double W;
    double s_min;
    double s_max;
    int left;
    int right;

    ParticleArray ps;        /* hat cuc bo + ghost */
    ParticleArray to_left;
    ParticleArray to_right;
    ParticleArray keep;
    ParticleArray gathered;
    IntArray neigh;
    Grid grid;

    int n_local;
    int step;
    int frame = 0;
    long long n_total = 0;
    long long n_local_ll;

    double t_comm = 0.0;
    double t_comp = 0.0;
    double t_wall0;
    double t_wall;
    double c0;
    double k0;

    double max_comm;
    double max_comp;
    double *all_comp = NULL;
    double *all_comm = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    params_default(&par);
    if (argc >= 2) {
        par.steps = atoi(argv[1]);
    }
    if (argc >= 3) {
        par.dx = atof(argv[2]);
    }
    if (argc >= 4) {
        par.dt = atof(argv[3]);
    }
    if (argc >= 5 && strcmp(argv[4], "pool") == 0) {
        params_set_pool(&par);
        use_pool = 1;
    }
    params_finalize(&par);

    W     = par.Lx / P;
    s_min = rank * W;
    s_max = (rank + 1) * W;
    left  = (rank == 0)     ? MPI_PROC_NULL : rank - 1;
    right = (rank == P - 1) ? MPI_PROC_NULL : rank + 1;

    pa_init(&ps);
    pa_init(&to_left);
    pa_init(&to_right);
    pa_init(&keep);
    pa_init(&gathered);
    ia_init(&neigh);

    /* Khoi tao: moi rank dung full khoi nuoc roi GIU hat thuoc lat cua minh. */
    {
        ParticleArray full;
        int i;
        pa_init(&full);
        sph_init_block(&full, &par);
        for (i = 0; i < full.size; i++) {
            int owner = (int) floor(full.data[i].pos.x / W);
            if (owner < 0) {
                owner = 0;
            }
            if (owner >= P) {
                owner = P - 1;
            }
            if (owner == rank) {
                pa_push(&ps, full.data[i]);
            }
        }
        pa_free(&full);
    }
    n_local = ps.size;

    n_local_ll = n_local;
    MPI_Reduce(&n_local_ll, &n_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("[mpi] P=%d ranks | N_tong=%lld hat | dx=%.3f h=%.3f | %d buoc | scenario=%s\n",
               P, n_total, par.dx, par.h, par.steps, use_pool ? "pool" : "dam");
    }

    grid_init(&grid, 1 << 16, par.h, vec3_zero());

    MPI_Barrier(MPI_COMM_WORLD);
    t_wall0 = MPI_Wtime();

    for (step = 0; step < par.steps; step++) {
        double wall_min;
        double wall_max;

        /* Bo ghost cua buoc truoc, chi giu hat cuc bo. */
        pa_resize(&ps, n_local);

        /* ---- (1) GHOST EXCHANGE lan 1: VI TRI (cho buoc tinh mat do) ---- */
        c0 = MPI_Wtime();
        build_halo(&ps, n_local, s_min, s_max, par.h, &to_left, &to_right);
        exchange(to_left.data,  to_left.size,  left,  right, &ps);
        exchange(to_right.data, to_right.size, right, left,  &ps);
        t_comm += MPI_Wtime() - c0;

        /* ---- (2) MAT DO + AP SUAT cho hat CUC BO ---- */
        k0 = MPI_Wtime();
        grid_build(&grid, ps.data, ps.size);
        sph_density_pressure(&ps, n_local, &grid, &par, &neigh);
        t_comp += MPI_Wtime() - k0;

        /* ---- GHOST EXCHANGE lan 2: rho,p (de hat ghost co ap suat dung) ---- */
        c0 = MPI_Wtime();
        pa_resize(&ps, n_local);
        build_halo(&ps, n_local, s_min, s_max, par.h, &to_left, &to_right);
        exchange(to_left.data,  to_left.size,  left,  right, &ps);
        exchange(to_right.data, to_right.size, right, left,  &ps);
        t_comm += MPI_Wtime() - c0;

        /* ---- (3) LUC + (4) TICH PHAN cho hat CUC BO ---- */
        k0 = MPI_Wtime();
        grid_build(&grid, ps.data, ps.size);
        sph_forces(&ps, n_local, &grid, &par, &neigh);
        wall_min = (rank == 0)     ? 0.0     : -1e30;
        wall_max = (rank == P - 1) ? par.Lx  :  1e30;
        sph_integrate(&ps, n_local, &par, wall_min, wall_max);
        t_comp += MPI_Wtime() - k0;

        /* ---- (5) MIGRATION: hat vuot bien lat -> rank lan can ---- */
        c0 = MPI_Wtime();
        {
            int i;
            pa_clear(&to_left);
            pa_clear(&to_right);
            pa_clear(&keep);
            for (i = 0; i < n_local; i++) {
                double x = ps.data[i].pos.x;
                if (x < s_min && left != MPI_PROC_NULL) {
                    pa_push(&to_left, ps.data[i]);
                } else if (x >= s_max && right != MPI_PROC_NULL) {
                    pa_push(&to_right, ps.data[i]);
                } else {
                    pa_push(&keep, ps.data[i]);
                }
            }
            pa_resize(&ps, 0);
            pa_append(&ps, keep.data, keep.size);
            exchange(to_left.data,  to_left.size,  left,  right, &ps);
            exchange(to_right.data, to_right.size, right, left,  &ps);
            n_local = ps.size;
        }
        t_comm += MPI_Wtime() - c0;

        /* ---- (6) OUTPUT: gather ve rank 0 ---- */
        if (step % par.out_every == 0) {
            int send_bytes = n_local * (int) sizeof(Particle);
            int *counts = NULL;
            int *displs = NULL;

            c0 = MPI_Wtime();
            if (rank == 0) {
                counts = (int *) malloc((size_t) P * sizeof(int));
                displs = (int *) malloc((size_t) P * sizeof(int));
            }
            MPI_Gather(&send_bytes, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                int total = 0;
                int i;
                for (i = 0; i < P; i++) {
                    displs[i] = total;
                    total += counts[i];
                }
                pa_resize(&gathered, total / (int) sizeof(Particle));
                MPI_Gatherv(ps.data, send_bytes, MPI_BYTE,
                            gathered.data, counts, displs, MPI_BYTE,
                            0, MPI_COMM_WORLD);
            } else {
                MPI_Gatherv(ps.data, send_bytes, MPI_BYTE,
                            NULL, NULL, NULL, MPI_BYTE, 0, MPI_COMM_WORLD);
            }
            t_comm += MPI_Wtime() - c0;

            if (rank == 0) {
                Stats st;
                char tag[32];
                io_write_frame("output", frame, gathered.data, gathered.size, W, P);
                io_write_vtk("output", frame, gathered.data, gathered.size);
                stats_compute(&st, gathered.data, gathered.size, params_mass(&par));
                snprintf(tag, sizeof(tag), "mpi f%03d", frame);
                stats_print(&st, tag);
                free(counts);
                free(displs);
            }
            frame++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t_wall = MPI_Wtime() - t_wall0;

    MPI_Reduce(&t_comm, &max_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_comp, &max_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        all_comp = (double *) malloc((size_t) P * sizeof(double));
        all_comm = (double *) malloc((size_t) P * sizeof(double));
    }
    MPI_Gather(&t_comp, 1, MPI_DOUBLE, all_comp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&t_comm, 1, MPI_DOUBLE, all_comm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double busiest = 0.0;
        double c_max = 0.0;
        double c_min = 1e30;
        int r;

        printf("[mpi] N_tong=%lld | P=%d\n", n_total, P);
        printf("[mpi] Tong thoi gian (co comm)      : %.3f s\n", t_wall);
        printf("[mpi] Compute (khong comm, max rank): %.3f s\n", max_comp);
        printf("[mpi] Communication (max rank)      : %.3f s (ghost+migration+gather)\n", max_comm);
        printf("[mpi] Comm/Total                    : %.1f%%\n", 100.0 * max_comm / t_wall);

        for (r = 0; r < P; r++) {
            double tot = all_comp[r] + all_comm[r];
            if (tot > busiest) {
                busiest = tot;
            }
        }
        printf("\n[per-rank] rank  compute(s)   comm(s)   total(s)   idle(s)\n");
        for (r = 0; r < P; r++) {
            double tot  = all_comp[r] + all_comm[r];
            double idle = busiest - tot;
            printf("           %4d  %9.3f  %8.3f  %8.3f  %7.3f\n",
                   r, all_comp[r], all_comm[r], tot, idle);
            if (all_comp[r] > c_max) {
                c_max = all_comp[r];
            }
            if (all_comp[r] < c_min) {
                c_min = all_comp[r];
            }
        }
        printf("[mpi] Mat can bang TAI (compute)    : %.1f%%  (=(max-min)/max)\n",
               c_max > 0.0 ? 100.0 * (c_max - c_min) / c_max : 0.0);

        io_write_timing("output/timing.csv", all_comp, all_comm, P);
        free(all_comp);
        free(all_comm);
    }

    grid_free(&grid);
    ia_free(&neigh);
    pa_free(&ps);
    pa_free(&to_left);
    pa_free(&to_right);
    pa_free(&keep);
    pa_free(&gathered);

    MPI_Finalize();
    return 0;
}
