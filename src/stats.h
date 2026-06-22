/* stats.h — Cac dai luong chan doan tren tap hat: dung de kiem tra dung dan
 * (trong tam, hop bao, mat do trung binh, toc do lon nhat, dong nang). */
#ifndef STATS_H
#define STATS_H

#include "particle.h"
#include "vec3.h"

typedef struct {
    int    n;
    Vec3   center;     /* trong tam */
    Vec3   bbox_min;   /* goc nho nhat hop bao */
    Vec3   bbox_max;   /* goc lon nhat hop bao */
    double rho_avg;    /* mat do trung binh */
    double rho_min;
    double rho_max;
    double speed_max;  /* |v| lon nhat */
    double kinetic;    /* tong dong nang (theo khoi luong moi hat) */
} Stats;

void stats_compute(Stats *s, const Particle *ps, int n, double mass);
void stats_print(const Stats *s, const char *tag);

#endif /* STATS_H */
