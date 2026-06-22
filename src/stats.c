/* stats.c — Cai dat tinh cac dai luong chan doan. */
#include "stats.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

void stats_compute(Stats *s, const Particle *ps, int n, double mass)
{
    int i;
    Vec3 sum = vec3_zero();

    s->n        = n;
    s->bbox_min = vec3_make(DBL_MAX, DBL_MAX, DBL_MAX);
    s->bbox_max = vec3_make(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    s->rho_avg  = 0.0;
    s->rho_min  = DBL_MAX;
    s->rho_max  = -DBL_MAX;
    s->speed_max = 0.0;
    s->kinetic   = 0.0;

    if (n <= 0) {
        s->center = vec3_zero();
        s->bbox_min = vec3_zero();
        s->bbox_max = vec3_zero();
        s->rho_min = 0.0;
        s->rho_max = 0.0;
        return;
    }

    for (i = 0; i < n; i++) {
        const Particle *a = &ps[i];
        double speed2 = vec3_norm2(a->vel);

        sum = vec3_add(sum, a->pos);

        if (a->pos.x < s->bbox_min.x) { s->bbox_min.x = a->pos.x; }
        if (a->pos.y < s->bbox_min.y) { s->bbox_min.y = a->pos.y; }
        if (a->pos.z < s->bbox_min.z) { s->bbox_min.z = a->pos.z; }
        if (a->pos.x > s->bbox_max.x) { s->bbox_max.x = a->pos.x; }
        if (a->pos.y > s->bbox_max.y) { s->bbox_max.y = a->pos.y; }
        if (a->pos.z > s->bbox_max.z) { s->bbox_max.z = a->pos.z; }

        s->rho_avg += a->rho;
        if (a->rho < s->rho_min) { s->rho_min = a->rho; }
        if (a->rho > s->rho_max) { s->rho_max = a->rho; }

        if (speed2 > s->speed_max) { s->speed_max = speed2; }  /* tam giu binh phuong */
        s->kinetic += 0.5 * mass * speed2;
    }

    s->center  = vec3_div(sum, (double) n);
    s->rho_avg /= (double) n;
    s->speed_max = sqrt(s->speed_max);   /* doi binh phuong -> toc do */
}

void stats_print(const Stats *s, const char *tag)
{
    printf("[%s] N=%d center=(%.4f,%.4f,%.4f) rho[avg=%.1f min=%.1f max=%.1f] vmax=%.3f KE=%.4f\n",
           tag, s->n,
           s->center.x, s->center.y, s->center.z,
           s->rho_avg, s->rho_min, s->rho_max,
           s->speed_max, s->kinetic);
}
