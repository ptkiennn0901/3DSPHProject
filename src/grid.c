/* grid.c — Cai dat bang bam khong gian. */
#include "grid.h"
#include <stdlib.h>
#include <math.h>

/* Toa do o cua mot diem theo tung truc. */
static int cell_coord(double v, double origin, double cell)
{
    return (int) floor((v - origin) / cell);
}

/* Bam (i,j,k) -> long long. Dich +512 de xu ly chi so am, goi 10 bit moi truc. */
static long long cell_hash(int i, int j, int k)
{
    long long ii = (long long) (i + 512) & 0x3FF;
    long long jj = (long long) (j + 512) & 0x3FF;
    long long kk = (long long) (k + 512) & 0x3FF;
    return (ii << 20) | (jj << 10) | kk;
}

/* Anh xa khoa long long -> chi so bucket. */
static int bucket_index(const Grid *g, long long key)
{
    unsigned long long h = (unsigned long long) key;
    h *= 1099511628211ULL;            /* tron bit kieu FNV */
    h ^= h >> 29;
    return (int) (h & (unsigned long long) (g->nbuckets - 1));
}

void grid_init(Grid *g, int nbuckets, double cell, Vec3 origin)
{
    /* Lam tron nbuckets len luy thua cua 2. */
    int n = 1;
    while (n < nbuckets) {
        n *= 2;
    }
    g->nbuckets = n;
    g->cell     = cell;
    g->origin   = origin;
    g->buckets  = (CellNode **) calloc((size_t) n, sizeof(CellNode *));
}

void grid_clear(Grid *g)
{
    int b;
    for (b = 0; b < g->nbuckets; b++) {
        CellNode *node = g->buckets[b];
        while (node != NULL) {
            CellNode *next = node->next;
            ia_free(&node->items);
            free(node);
            node = next;
        }
        g->buckets[b] = NULL;
    }
}

void grid_free(Grid *g)
{
    grid_clear(g);
    free(g->buckets);
    g->buckets  = NULL;
    g->nbuckets = 0;
}

/* Tim node cua khoa trong bucket, hoac tao moi neu chua co. */
static CellNode *grid_get_or_create(Grid *g, long long key)
{
    int b = bucket_index(g, key);
    CellNode *node = g->buckets[b];
    while (node != NULL) {
        if (node->key == key) {
            return node;
        }
        node = node->next;
    }
    node = (CellNode *) malloc(sizeof(CellNode));
    node->key = key;
    ia_init(&node->items);
    node->next = g->buckets[b];
    g->buckets[b] = node;
    return node;
}

/* Tim node theo khoa (chi doc), tra NULL neu khong co. */
static const CellNode *grid_find(const Grid *g, long long key)
{
    int b = bucket_index(g, key);
    const CellNode *node = g->buckets[b];
    while (node != NULL) {
        if (node->key == key) {
            return node;
        }
        node = node->next;
    }
    return NULL;
}

void grid_build(Grid *g, const Particle *ps, int n)
{
    int i;
    grid_clear(g);
    for (i = 0; i < n; i++) {
        int ci = cell_coord(ps[i].pos.x, g->origin.x, g->cell);
        int cj = cell_coord(ps[i].pos.y, g->origin.y, g->cell);
        int ck = cell_coord(ps[i].pos.z, g->origin.z, g->cell);
        long long key = cell_hash(ci, cj, ck);
        CellNode *node = grid_get_or_create(g, key);
        ia_push(&node->items, i);
    }
}

void grid_neighbors(const Grid *g, Vec3 p, IntArray *out)
{
    int ci = cell_coord(p.x, g->origin.x, g->cell);
    int cj = cell_coord(p.y, g->origin.y, g->cell);
    int ck = cell_coord(p.z, g->origin.z, g->cell);
    int di, dj, dk;
    ia_clear(out);
    for (di = -1; di <= 1; di++) {
        for (dj = -1; dj <= 1; dj++) {
            for (dk = -1; dk <= 1; dk++) {
                long long key = cell_hash(ci + di, cj + dj, ck + dk);
                const CellNode *node = grid_find(g, key);
                if (node != NULL) {
                    ia_append(out, node->items.data, node->items.size);
                }
            }
        }
    }
}
