# Makefile — build 3 chuong trinh C: hello (test MPI), sph_serial, sph_mpi
CC       = cc
MPICC    = mpicc
CFLAGS   = -O2 -std=gnu11 -Wall -Wextra
LDFLAGS  = -lm

CORE_SRC = src/vec3.c src/darray.c src/kernels.c src/params.c src/grid.c \
           src/sph_core.c src/io.c src/stats.c

all: hello sph_serial sph_mpi

hello: src/hello_mpi.c
	$(MPICC) $(CFLAGS) -o hello src/hello_mpi.c

sph_serial: src/sph_serial.c $(CORE_SRC)
	$(CC) $(CFLAGS) -o sph_serial src/sph_serial.c $(CORE_SRC) $(LDFLAGS)

sph_mpi: src/sph_mpi.c $(CORE_SRC)
	$(MPICC) $(CFLAGS) -o sph_mpi src/sph_mpi.c $(CORE_SRC) $(LDFLAGS)

clean:
	rm -f hello sph_serial sph_mpi output/frame_*.csv output/timing.csv

.PHONY: all clean
