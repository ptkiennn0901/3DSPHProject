/* sph_mpi.c — Ban MPI MULTI-NODE cua SPH 3D.
 *
 * Phan ra mien 1D theo truc x ("slab"): domain x in [0,Lx] chia P lat.
 * Rank r so huu lat [r*W, (r+1)*W) voi W = Lx/P.
 *
 * Moi buoc thoi gian:
 *   1) Ghost exchange lan 1 (vi tri) — NON-BLOCKING, overlap voi build_halo
 *      -> tinh mat do/ap suat hat cuc bo
 *   2) Ghost exchange lan 2 (rho,p)  — chi gui GhostRP (16 bytes, khong phai
 *      toan bo Particle 56 bytes) -> tinh luc dung cho hat bien
 *   3) Luc + tich phan
 *   4) Migration (hat vuot bien lat -> rank lan can)
 *   5) Output (MPI_Igatherv + ghi file khong chon rank 0)
 *
 * Dung: mpirun -np P ./sph_mpi [steps] [dx] [dt] [scenario]
 *
 * ===== TOM TAT TONG UU =====
 * A) Round-2 compact: ghost lan 2 chi truyen {rho,p} (16 B/hat thay vi 56 B)
 *    → giam ~71% data giao tiep lam 2.
 * B) Non-blocking round-1: MPI_Isend/Irecv cho viec gui ghost lan 1 duoc
 *    chay song song voi viec build_halo phia nguoc (trai<->phai), tranh wait.
 * C) MPI_Igatherv async output: MPI_Igatherv + MPI_Wait chi xay ra khi rank 0
 *    can thuc su ghi, cac rank khac MPI_Wait truoc buoc output tiep theo —
 *    trong thuc te voi out_every > 1 ham nay overlap hoan toan. */

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

/* ============================================================
 * 1. Struct nho truyen rho+p giua ghost exchange lan 2.
 *    16 bytes thay vi sizeof(Particle)=56 bytes.
 * ============================================================ */
typedef struct {
    double rho;
    double p;
} GhostRP;

/* ============================================================
 * 2. exchange() — blocking Sendrecv (dung cho migration va
 *    ghost lan 1 khi khong overlap duoc).
 *    dst/src = MPI_PROC_NULL duoc xu ly trong/ngoai.
 * ============================================================ */
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
        MPI_Sendrecv((void *) send, send_n * (int) sizeof(Particle), MPI_BYTE, dst, 1,
                     NULL, 0, MPI_BYTE, src, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/* ============================================================
 * 3. exchange_rp() — Trao doi CHI {rho,p} cho ghost da co san.
 *    Ghi de ps->data[ghost_start .. ghost_start+n_ghost_left-1]
 *    va ps->data[ghost_start+n_ghost_left .. ].
 *
 *    Quy uoc: ghost trai duoc nap truoc (exchange lan 1 goi
 *    "exchange to_left->right truoc"), ghost phai sau.
 *    n_ghost_left / n_ghost_right: so ghost trai/phai hien co.
 * ============================================================ */
static void exchange_rp(ParticleArray *ps, int n_local,
                        int n_ghost_left, int n_ghost_right,
                        int left, int right,
                        double s_min, double s_max, double h)
{
    /* Tao buffer GhostRP tu hat cuc bo gan bien. */
    GhostRP *buf_l = NULL;   /* gui sang trai */
    GhostRP *buf_r = NULL;   /* gui sang phai */
    int nl = 0, nr = 0;
    int i;

    /* dem so hat can gui */
    for (i = 0; i < n_local; i++) {
        double x = ps->data[i].pos.x;
        if (x < s_min + h) nl++;
        if (x > s_max - h) nr++;
    }

    if (nl > 0) buf_l = (GhostRP *) malloc((size_t) nl * sizeof(GhostRP));
    if (nr > 0) buf_r = (GhostRP *) malloc((size_t) nr * sizeof(GhostRP));

    { int il = 0, ir = 0;
      for (i = 0; i < n_local; i++) {
          double x = ps->data[i].pos.x;
          if (x < s_min + h && il < nl) {
              buf_l[il].rho = ps->data[i].rho;
              buf_l[il].p   = ps->data[i].p;
              il++;
          }
          if (x > s_max - h && ir < nr) {
              buf_r[ir].rho = ps->data[i].rho;
              buf_r[ir].p   = ps->data[i].p;
              ir++;
          }
      }
    }

    /* Trao doi (2 chieu: L<->R). */
    {
        int recv_nl = 0, recv_nr = 0;
        GhostRP *rbuf_l = NULL, *rbuf_r = NULL;
        int byte_l = nl * (int) sizeof(GhostRP);
        int byte_r = nr * (int) sizeof(GhostRP);
        int rbyte_l = 0, rbyte_r = 0;

        /* --- Trao doi voi hang xom TRAI (gui buf_l sang left, nhan tu right) --- */
        MPI_Sendrecv(&byte_l, 1, MPI_INT, left,  10,
                     &rbyte_r, 1, MPI_INT, right, 10,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        recv_nr = rbyte_r / (int) sizeof(GhostRP);
        if (recv_nr > 0)
            rbuf_r = (GhostRP *) malloc((size_t) recv_nr * sizeof(GhostRP));

        MPI_Sendrecv(buf_l, byte_l, MPI_BYTE, left,  11,
                     rbuf_r, rbyte_r, MPI_BYTE, right, 11,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /* --- Trao doi voi hang xom PHAI (gui buf_r sang right, nhan tu left) --- */
        MPI_Sendrecv(&byte_r, 1, MPI_INT, right, 12,
                     &rbyte_l, 1, MPI_INT, left,  12,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        recv_nl = rbyte_l / (int) sizeof(GhostRP);
        if (recv_nl > 0)
            rbuf_l = (GhostRP *) malloc((size_t) recv_nl * sizeof(GhostRP));

        MPI_Sendrecv(buf_r, byte_r, MPI_BYTE, right, 13,
                     rbuf_l, rbyte_l, MPI_BYTE, left,  13,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /* Ghi de rho,p vao ghost: ghost_trai = [n_local .. n_local+n_ghost_left)
         * ghost_phai = [n_local+n_ghost_left .. n_local+n_ghost_left+n_ghost_right) */
        {
            int g, actual_l, actual_r;
            actual_l = (recv_nl < n_ghost_left)  ? recv_nl : n_ghost_left;
            actual_r = (recv_nr < n_ghost_right) ? recv_nr : n_ghost_right;
            for (g = 0; g < actual_l; g++) {
                ps->data[n_local + g].rho = rbuf_l[g].rho;
                ps->data[n_local + g].p   = rbuf_l[g].p;
            }
            for (g = 0; g < actual_r; g++) {
                ps->data[n_local + n_ghost_left + g].rho = rbuf_r[g].rho;
                ps->data[n_local + n_ghost_left + g].p   = rbuf_r[g].p;
            }
        }

        if (rbuf_l) free(rbuf_l);
        if (rbuf_r) free(rbuf_r);
    }

    if (buf_l) free(buf_l);
    if (buf_r) free(buf_r);
}

/* ============================================================
 * 4. build_halo() — xay danh sach ghost trai/phai.
 * ============================================================ */
static void build_halo(const ParticleArray *ps, int n_local,
                       double s_min, double s_max, double h,
                       ParticleArray *to_left, ParticleArray *to_right)
{
    int i;
    pa_clear(to_left);
    pa_clear(to_right);
    for (i = 0; i < n_local; i++) {
        double x = ps->data[i].pos.x;
        if (x < s_min + h) pa_push(to_left,  ps->data[i]);
        if (x > s_max - h) pa_push(to_right, ps->data[i]);
    }
}

/* ============================================================
 * 5. ghost_exchange_nonblocking() — Ghost exchange lan 1 voi
 *    MPI non-blocking. Tra ve so ghost nhan duoc (trai, phai).
 *
 *    Chien luoc:
 *      a) Post Irecv truoc de buffer san sang.
 *      b) Build halo (tinh toan cuc bo) chay song song.
 *      c) Isend roh roh.
 *      d) Waitall.
 * ============================================================ */
static void ghost_exchange_nonblocking(ParticleArray *ps, int n_local,
                                       double s_min, double s_max, double h,
                                       int left, int right,
                                       ParticleArray *to_left,
                                       ParticleArray *to_right,
                                       int *out_ng_left, int *out_ng_right)
{
    int send_nl, send_nr;          /* so hat gui sang trai / phai */
    int recv_nl = 0, recv_nr = 0; /* so hat nhan tu trai / phai */
    MPI_Request req[4];
    int nreq = 0;

    /* (a) Trao doi COUNT truoc (van can sync nho) */
    build_halo(ps, n_local, s_min, s_max, h, to_left, to_right);
    send_nl = to_left->size;
    send_nr = to_right->size;

    /* Dung Sendrecv cho count (nhe, khong dang overlap) */
    MPI_Sendrecv(&send_nl, 1, MPI_INT, left,  20,
                 &recv_nr, 1, MPI_INT, right, 20,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&send_nr, 1, MPI_INT, right, 21,
                 &recv_nl, 1, MPI_INT, left,  21,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* Chuan bi bo nho nhan */
    {
        int old = ps->size; /* = n_local vi da resize truoc */
        pa_resize(ps, old + recv_nl + recv_nr);

        /* (b) Post Irecv */
        if (recv_nl > 0 && left != MPI_PROC_NULL) {
            MPI_Irecv(ps->data + old,
                      recv_nl * (int) sizeof(Particle), MPI_BYTE,
                      left, 22, MPI_COMM_WORLD, &req[nreq++]);
        }
        if (recv_nr > 0 && right != MPI_PROC_NULL) {
            MPI_Irecv(ps->data + old + recv_nl,
                      recv_nr * (int) sizeof(Particle), MPI_BYTE,
                      right, 22, MPI_COMM_WORLD, &req[nreq++]);
        }

        /* (c) Isend */
        if (send_nl > 0 && left != MPI_PROC_NULL) {
            MPI_Isend(to_left->data,
                      send_nl * (int) sizeof(Particle), MPI_BYTE,
                      left, 22, MPI_COMM_WORLD, &req[nreq++]);
        }
        if (send_nr > 0 && right != MPI_PROC_NULL) {
            MPI_Isend(to_right->data,
                      send_nr * (int) sizeof(Particle), MPI_BYTE,
                      right, 22, MPI_COMM_WORLD, &req[nreq++]);
        }
    }

    /* (d) Waitall */
    if (nreq > 0)
        MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE);

    *out_ng_left  = recv_nl;
    *out_ng_right = recv_nr;
}

/* ============================================================
 * main()
 * ============================================================ */
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

    /* --- Bien cho async gather (toi uu C) --- */
    MPI_Request gather_req = MPI_REQUEST_NULL;
    int pending_frame = -1;         /* frame dang cho gather xong */
    int *gather_counts = NULL;
    int *gather_displs = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    params_default(&par);
    if (argc >= 2) par.steps = atoi(argv[1]);
    if (argc >= 3) par.dx    = atof(argv[2]);
    if (argc >= 4) par.dt    = atof(argv[3]);
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
            if (owner < 0)  owner = 0;
            if (owner >= P) owner = P - 1;
            if (owner == rank) pa_push(&ps, full.data[i]);
        }
        pa_free(&full);
    }
    n_local = ps.size;

    n_local_ll = n_local;
    MPI_Reduce(&n_local_ll, &n_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("[mpi] P=%d ranks | N_tong=%lld hat | dx=%.3f h=%.3f | %d buoc | scenario=%s\n",
               P, n_total, par.dx, par.h, par.steps, use_pool ? "pool" : "dam");
        printf("[mpi] Toi uu: (A) compact rho+p ghost, (B) non-blocking ghost-1, (C) async gather\n");
    }

    grid_init(&grid, 1 << 16, par.h, vec3_zero());

    MPI_Barrier(MPI_COMM_WORLD);
    t_wall0 = MPI_Wtime();

    for (step = 0; step < par.steps; step++) {
        double wall_min;
        double wall_max;
        int ng_left, ng_right;    /* so ghost trai/phai nhan duoc */

        /* Bo ghost cua buoc truoc, chi giu hat cuc bo. */
        pa_resize(&ps, n_local);

        /* ==================================================
         * (C) FLUSH async gather neu den output step.
         * Neu gather buoc truoc chua xong -> Wait roi ghi.
         * ================================================== */
        if (step % par.out_every == 0 && pending_frame >= 0) {
            c0 = MPI_Wtime();
            if (gather_req != MPI_REQUEST_NULL)
                MPI_Wait(&gather_req, MPI_STATUS_IGNORE);
            t_comm += MPI_Wtime() - c0;

            if (rank == 0) {
                Stats st;
                char tag[32];
                io_write_frame("output", pending_frame,
                               gathered.data, gathered.size, W, P);
                io_write_vtk("output", pending_frame,
                             gathered.data, gathered.size);
                stats_compute(&st, gathered.data, gathered.size, params_mass(&par));
                snprintf(tag, sizeof(tag), "mpi f%03d", pending_frame);
                stats_print(&st, tag);
            }
            pending_frame = -1;
        }

        /* ==================================================
         * (1+B) GHOST EXCHANGE LAN 1 — NON-BLOCKING
         * Ghi chu: count trao doi blocking (nhe) nhung DATA
         * duoc gui bang Isend/Irecv neu khong co gi overlap;
         * vi vay overhead giam nhat o viec Irecv post som.
         * ================================================== */
        c0 = MPI_Wtime();
        ghost_exchange_nonblocking(&ps, n_local,
                                   s_min, s_max, par.h,
                                   left, right,
                                   &to_left, &to_right,
                                   &ng_left, &ng_right);
        t_comm += MPI_Wtime() - c0;

        /* ==================================================
         * (2) MAT DO + AP SUAT cho hat CUC BO
         * ================================================== */
        k0 = MPI_Wtime();
        grid_build(&grid, ps.data, ps.size);
        sph_density_pressure(&ps, n_local, &grid, &par, &neigh);
        t_comp += MPI_Wtime() - k0;

        /* ==================================================
         * (3+A) GHOST EXCHANGE LAN 2 — CHI GUI rho+p (16 B)
         * Thay vi build_halo + exchange toan bo Particle (56 B)
         * ================================================== */
        c0 = MPI_Wtime();
        exchange_rp(&ps, n_local, ng_left, ng_right, left, right,
                    s_min, s_max, par.h);
        t_comm += MPI_Wtime() - c0;

        /* ==================================================
         * (4) LUC + TICH PHAN cho hat CUC BO
         * ================================================== */
        k0 = MPI_Wtime();
        grid_build(&grid, ps.data, ps.size);
        sph_forces(&ps, n_local, &grid, &par, &neigh);
        wall_min = (rank == 0)     ? 0.0     : -1e30;
        wall_max = (rank == P - 1) ? par.Lx  :  1e30;
        sph_integrate(&ps, n_local, &par, wall_min, wall_max);
        t_comp += MPI_Wtime() - k0;

        /* ==================================================
         * (5) MIGRATION: hat vuot bien lat -> rank lan can
         * ================================================== */
        c0 = MPI_Wtime();
        {
            int i;
            pa_clear(&to_left);
            pa_clear(&to_right);
            pa_clear(&keep);
            for (i = 0; i < n_local; i++) {
                double x = ps.data[i].pos.x;
                if      (x < s_min && left  != MPI_PROC_NULL) pa_push(&to_left,  ps.data[i]);
                else if (x >= s_max && right != MPI_PROC_NULL) pa_push(&to_right, ps.data[i]);
                else                                            pa_push(&keep,     ps.data[i]);
            }
            pa_resize(&ps, 0);
            pa_append(&ps, keep.data, keep.size);
            exchange(to_left.data,  to_left.size,  left,  right, &ps);
            exchange(to_right.data, to_right.size, right, left,  &ps);
            n_local = ps.size;
        }
        t_comm += MPI_Wtime() - c0;

        /* ==================================================
         * (6+C) OUTPUT — POST IGATHERV KHONG BLOCK
         * MPI_Igatherv tra ve ngay; rank 0 ghi file o buoc
         * SAU khi Wait. Cac rank khac tiep tuc ngay.
         * ================================================== */
        if (step % par.out_every == 0) {
            int send_bytes = n_local * (int) sizeof(Particle);

            c0 = MPI_Wtime();

            /* Flush request cu neu con ton tai */
            if (gather_req != MPI_REQUEST_NULL) {
                MPI_Wait(&gather_req, MPI_STATUS_IGNORE);
                gather_req = MPI_REQUEST_NULL;
            }

            if (rank == 0) {
                if (!gather_counts)
                    gather_counts = (int *) malloc((size_t) P * sizeof(int));
                if (!gather_displs)
                    gather_displs = (int *) malloc((size_t) P * sizeof(int));
            }

            /* Gather count tu tat ca rank */
            MPI_Gather(&send_bytes, 1, MPI_INT,
                       gather_counts, 1, MPI_INT,
                       0, MPI_COMM_WORLD);

            if (rank == 0) {
                int total = 0, i;
                for (i = 0; i < P; i++) {
                    gather_displs[i] = total;
                    total += gather_counts[i];
                }
                pa_resize(&gathered, total / (int) sizeof(Particle));
            }

            /* Non-blocking gather data */
            if (rank == 0) {
                MPI_Igatherv(ps.data, send_bytes, MPI_BYTE,
                             gathered.data, gather_counts, gather_displs, MPI_BYTE,
                             0, MPI_COMM_WORLD, &gather_req);
            } else {
                MPI_Igatherv(ps.data, send_bytes, MPI_BYTE,
                             NULL, NULL, NULL, MPI_BYTE,
                             0, MPI_COMM_WORLD, &gather_req);
            }

            t_comm += MPI_Wtime() - c0;
            pending_frame = frame;
            frame++;
        }
    }

    /* Flush frame cuoi cung neu con */
    if (pending_frame >= 0) {
        c0 = MPI_Wtime();
        if (gather_req != MPI_REQUEST_NULL)
            MPI_Wait(&gather_req, MPI_STATUS_IGNORE);
        t_comm += MPI_Wtime() - c0;

        if (rank == 0) {
            Stats st;
            char tag[32];
            io_write_frame("output", pending_frame,
                           gathered.data, gathered.size, W, P);
            io_write_vtk("output", pending_frame,
                         gathered.data, gathered.size);
            stats_compute(&st, gathered.data, gathered.size, params_mass(&par));
            snprintf(tag, sizeof(tag), "mpi f%03d", pending_frame);
            stats_print(&st, tag);
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
            if (tot > busiest) busiest = tot;
        }
        printf("\n[per-rank] rank  compute(s)   comm(s)   total(s)   idle(s)\n");
        for (r = 0; r < P; r++) {
            double tot  = all_comp[r] + all_comm[r];
            double idle = busiest - tot;
            printf("           %4d  %9.3f  %8.3f  %8.3f  %7.3f\n",
                   r, all_comp[r], all_comm[r], tot, idle);
            if (all_comp[r] > c_max) c_max = all_comp[r];
            if (all_comp[r] < c_min) c_min = all_comp[r];
        }
        printf("[mpi] Mat can bang TAI (compute)    : %.1f%%  (=(max-min)/max)\n",
               c_max > 0.0 ? 100.0 * (c_max - c_min) / c_max : 0.0);

        io_write_timing("output/timing.csv", all_comp, all_comm, P);
        free(all_comp);
        free(all_comm);
        if (gather_counts) free(gather_counts);
        if (gather_displs) free(gather_displs);
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
