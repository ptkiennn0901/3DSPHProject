/* io.h — Ghi frame ra CSV de Python dung hoat anh 3D, va ghi timing. */
#ifndef IO_H
#define IO_H

#include "particle.h"

/* Ghi n hat dau cua mang ps ra output/frame_XXXX.csv.
 * own_axis_width > 0: tinh cot 'rank' theo so huu lat x (de to mau). */
void io_write_frame(const char *dir, int frame,
                    const Particle *ps, int n, double own_axis_width, int nproc);

/* Ghi bang thoi gian tung rank ra output/timing.csv (cho TN2). */
void io_write_timing(const char *path, const double *compute,
                     const double *comm, int nproc);

/* Ghi frame dang VTK legacy (POLYDATA) de mo bang ParaView. */
void io_write_vtk(const char *dir, int frame, const Particle *ps, int n);

#endif /* IO_H */
