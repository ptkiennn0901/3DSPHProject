/* darray.c — Cai dat mang dong. Tang dung luong gap doi khi day. */
#include "darray.h"
#include <stdlib.h>
#include <string.h>

/* ================= IntArray ================= */

void ia_init(IntArray *a)
{
    a->data = NULL;
    a->size = 0;
    a->cap  = 0;
}

void ia_reserve(IntArray *a, int cap)
{
    if (cap <= a->cap) {
        return;
    }
    int newcap = a->cap > 0 ? a->cap : 8;
    while (newcap < cap) {
        newcap *= 2;
    }
    int *p = (int *) realloc(a->data, (size_t) newcap * sizeof(int));
    if (p == NULL) {
        return;
    }
    a->data = p;
    a->cap  = newcap;
}

void ia_push(IntArray *a, int value)
{
    if (a->size >= a->cap) {
        ia_reserve(a, a->size + 1);
    }
    a->data[a->size] = value;
    a->size += 1;
}

void ia_clear(IntArray *a)
{
    a->size = 0;
}

void ia_free(IntArray *a)
{
    free(a->data);
    a->data = NULL;
    a->size = 0;
    a->cap  = 0;
}

void ia_append(IntArray *a, const int *src, int n)
{
    if (n <= 0) {
        return;
    }
    ia_reserve(a, a->size + n);
    memcpy(a->data + a->size, src, (size_t) n * sizeof(int));
    a->size += n;
}

/* ================= ParticleArray ================= */

void pa_init(ParticleArray *a)
{
    a->data = NULL;
    a->size = 0;
    a->cap  = 0;
}

void pa_reserve(ParticleArray *a, int cap)
{
    if (cap <= a->cap) {
        return;
    }
    int newcap = a->cap > 0 ? a->cap : 8;
    while (newcap < cap) {
        newcap *= 2;
    }
    Particle *p = (Particle *) realloc(a->data, (size_t) newcap * sizeof(Particle));
    if (p == NULL) {
        return;
    }
    a->data = p;
    a->cap  = newcap;
}

void pa_push(ParticleArray *a, Particle value)
{
    if (a->size >= a->cap) {
        pa_reserve(a, a->size + 1);
    }
    a->data[a->size] = value;
    a->size += 1;
}

void pa_clear(ParticleArray *a)
{
    a->size = 0;
}

void pa_resize(ParticleArray *a, int n)
{
    if (n < 0) {
        n = 0;
    }
    if (n > a->size) {
        pa_reserve(a, n);
    }
    a->size = n;
}

void pa_free(ParticleArray *a)
{
    free(a->data);
    a->data = NULL;
    a->size = 0;
    a->cap  = 0;
}

void pa_append(ParticleArray *a, const Particle *src, int n)
{
    if (n <= 0) {
        return;
    }
    pa_reserve(a, a->size + n);
    memcpy(a->data + a->size, src, (size_t) n * sizeof(Particle));
    a->size += n;
}
