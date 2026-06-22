/* darray.h — Mang dong tu cai (thay cho std::vector cua C++).
 * IntArray      : danh sach chi so hat (cho lan can va o luoi).
 * ParticleArray : danh sach hat (local + ghost). */
#ifndef DARRAY_H
#define DARRAY_H

#include "particle.h"

/* ---- Mang dong cac so nguyen ---- */
typedef struct {
    int *data;
    int  size;
    int  cap;
} IntArray;

void ia_init(IntArray *a);
void ia_reserve(IntArray *a, int cap);
void ia_push(IntArray *a, int value);
void ia_clear(IntArray *a);          /* size = 0, giu bo nho */
void ia_free(IntArray *a);
void ia_append(IntArray *a, const int *src, int n);

/* ---- Mang dong cac hat ---- */
typedef struct {
    Particle *data;
    int       size;
    int       cap;
} ParticleArray;

void pa_init(ParticleArray *a);
void pa_reserve(ParticleArray *a, int cap);
void pa_push(ParticleArray *a, Particle value);
void pa_clear(ParticleArray *a);
void pa_resize(ParticleArray *a, int n);   /* dat size = n (n <= size) */
void pa_free(ParticleArray *a);
void pa_append(ParticleArray *a, const Particle *src, int n);

#endif /* DARRAY_H */
