/* grid.h — Luoi khong gian deu de tim lan can voi do phuc tap O(N).
 * O lap phuong canh h; moi hat chi quet 27 o (3x3x3) xung quanh.
 *
 * Vi C khong co unordered_map, ta tu cai bang bam (hash table) dung
 * chaining: moi bucket la mot danh sach lien ket cac CellNode, moi
 * CellNode giu khoa o (i,j,k) bam thanh long long va danh sach chi so hat. */
#ifndef GRID_H
#define GRID_H

#include "vec3.h"
#include "darray.h"

typedef struct CellNode {
    long long        key;     /* khoa o (i,j,k) da bam */
    IntArray         items;   /* chi so cac hat trong o */
    struct CellNode *next;    /* nut tiep theo cung bucket */
} CellNode;

typedef struct {
    CellNode **buckets;       /* mang con tro dau bucket */
    int        nbuckets;      /* so bucket (luy thua cua 2) */
    double     cell;          /* canh o = h */
    Vec3       origin;        /* goc nho nhat cua luoi */
} Grid;

void grid_init(Grid *g, int nbuckets, double cell, Vec3 origin);
void grid_clear(Grid *g);                       /* giai phong cac node, giu mang bucket */
void grid_free(Grid *g);
void grid_build(Grid *g, const Particle *ps, int n);
void grid_neighbors(const Grid *g, Vec3 p, IntArray *out);

#endif /* GRID_H */
