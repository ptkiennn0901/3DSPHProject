/* sph_core.h — Loi vat ly SPH 3D. Bon buoc:
 *   1) Mat do:    rho_i = sum_j m_j W_poly6(|r_i-r_j|, h)
 *   2) Ap suat:   p_i   = k (rho_i - rho0)
 *   3) Luc:       F = F_ap_suat + F_nhot + trong_luc  (mat do luc)
 *   4) Tich phan: a = F/rho; v += dt a; r += dt v + phan xa 6 tuong
 *
 * Quy uoc MPI: mang hat chua hat CUC BO o [0,n_local) roi toi hat GHOST.
 * Cap nhat hat cuc bo nhung quet lan can tren toan mang (ke ca ghost). */
#ifndef SPH_CORE_H
#define SPH_CORE_H

#include "params.h"
#include "darray.h"
#include "grid.h"

/* Tao khoi nuoc (dam o goc hoac pool day be) vao mang 'out'. */
void sph_init_block(ParticleArray *out, const Params *P);

/* Buoc 1+2: mat do va ap suat cho n_local hat dau. */
void sph_density_pressure(ParticleArray *ps, int n_local,
                          const Grid *grid, const Params *P, IntArray *neigh);

/* Buoc 3: luc ap suat + nhot + trong luc cho n_local hat dau. */
void sph_forces(ParticleArray *ps, int n_local,
                const Grid *grid, const Params *P, IntArray *neigh);

/* Buoc 4: tich phan + phan xa tuong cho n_local hat dau.
 * wall_xmin/wall_xmax: tuong that theo truc x (serial: 0..Lx). */
void sph_integrate(ParticleArray *ps, int n_local, const Params *P,
                   double wall_xmin, double wall_xmax);

/* Mat do tham chieu O(N^2) (khong dung luoi) — de kiem tra cheo do dung. */
double sph_density_bruteforce(const ParticleArray *ps, int i, const Params *P);

/* Tu kiem tra: so sanh mat do tinh bang luoi vs O(N^2) tren vai hat mau.
 * Tra ve sai so tuong doi lon nhat (cang nho cang tot). */
double sph_self_test(ParticleArray *ps, const Grid *grid,
                     const Params *P, IntArray *neigh, int samples);

#endif /* SPH_CORE_H */
