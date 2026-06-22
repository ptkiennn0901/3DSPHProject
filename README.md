# 3D SPH Fluid Simulation — Parallel (OpenMPI, multi-machine) — bản C

Mô phỏng chất lưu SPH 3D (dam break) song song bằng **MPI**, viết hoàn toàn bằng **C** (C11/gnu11).
Có bản **serial** (nền tảng) và bản **MPI multi-node** dùng chung lõi vật lý. Visualize bằng Python.

Vì C không có `std::vector` / `unordered_map`, dự án **tự cài** mảng động (`darray`) và bảng băm
không gian (`grid`) — đây cũng là phần code đáng kể của bản C.

## Cấu trúc (src/)

| File | Vai trò |
|------|---------|
| `vec3.h/.c`      | toán vector 3D (add, sub, dot, cross, norm, lerp, ...) |
| `particle.h`     | struct hạt (POD, gửi qua MPI bằng MPI_BYTE) |
| `darray.h/.c`    | mảng động IntArray / ParticleArray (thay std::vector) |
| `grid.h/.c`      | bảng băm không gian, tìm lân cận O(N) qua 27 ô |
| `kernels.h/.c`   | kernel SPH 3D: Poly6 / Spiky grad / Viscosity lap |
| `params.h/.c`    | tham số mô phỏng + kịch bản dam/pool |
| `sph_core.h/.c`  | 4 bước SPH + self-test lưới vs O(N²) |
| `stats.h/.c`     | chẩn đoán: trọng tâm, hộp bao, mật độ, dộng năng |
| `io.h/.c`        | ghi frame CSV, timing CSV, và **VTK cho ParaView** |
| `sph_serial.c`   | bản 1 tiến trình |
| `sph_mpi.c`      | bản MPI: domain decomp + ghost 2 pha + migration + gather |
| `hello_mpi.c`    | test cluster |

## Build

```bash
make            # cc + mpicc, tạo: hello, sph_serial, sph_mpi
```

## Chạy thử nhanh trên 1 máy

```bash
./sph_serial 300 0.04 0.0004                 # ~700 hạt, kịch bản dam
mpirun --oversubscribe -np 3 ./sph_mpi 300 0.04 0.0004
```
Bản serial in **self-test** (mật độ lưới vs O(N²), phải PASS ~1e-16), và mỗi frame in `stats`
(trọng tâm, mật độ, vmax, động năng). Trọng tâm của serial và MPI phải trùng nhau.

Tham số: `./sph_xxx [steps] [dx] [dt] [scenario]`
- `steps` số bước; `dx` khoảng cách hạt (N ≈ thể_tích/dx³); `dt` bước thời gian (≈ 0.01·dx);
- `scenario` = `dam` (mặc định, khối góc) hoặc `pool` (đầy bể, dùng benchmark).

## Tăng số hạt

N do `dx` quyết định. Giảm `dx` một nửa → N tăng 8 lần. **Phải giảm dt theo (dt ≈ 0.01·dx)**
nếu không mô phỏng mất ổn định (mật độ tụt khỏi 1000).

| dx | N (dam) | N (pool) | dt |
|----|---------|----------|-----|
| 0.02  | 6K   | 32K  | 0.0002  |
| 0.015 | 13.5K| 75K  | 0.00015 |
| 0.012 | 27K  | 146K | 0.00012 |
| 0.010 | 48K  | 261K | 0.0001  |

## Chạy trên cụm 3 máy thật

Cài trên cả 3 máy: `build-essential openmpi-bin libopenmpi-dev` (cùng phiên bản). Cấu hình SSH
không mật khẩu, `/etc/hosts`, đồng bộ binary (NFS/rsync, cùng đường dẫn), tạo `~/hostfile`
(xem `hostfile.example`). Test: `mpirun --hostfile ~/hostfile -np 6 ./hello`. Rồi:
```bash
mpirun --hostfile ~/hostfile -np 6 ./sph_mpi 2000 0.02 0.0002 pool
```

## Trực quan hóa

```bash
pip3 install numpy matplotlib
python3 scripts/visualize.py --color rank        # GIF tô màu theo rank/máy
python3 scripts/visualize.py --frame 20 --color speed
```
Hoặc mở các file `output/frame_*.vtk` bằng **ParaView** (chất lượng cao hơn cho N lớn).

## Bộ thí nghiệm tự động (TN0–TN3 cho báo cáo)

```bash
export CORES=12 HOSTFILE=~/hostfile STEPS=1000
./scripts/run_experiments.sh tn0 0.02 0.0002      # đúng đắn serial vs MPI
./scripts/run_experiments.sh tn1                  # quét N -> tìm N cho 2-3 phút
./scripts/run_experiments.sh tn2 <dx*> <dt*>      # cân bằng tải
./scripts/run_experiments.sh tn3 <dx_2N> <dt>     # speedup
python3 scripts/plot_results.py                   # vẽ results/*.png
```
Mọi benchmark dùng kịch bản **pool** (cân bằng tải, đã đặt sẵn trong script).

## Thiết kế song song (tóm tắt cho báo cáo)

- **Cấp độ:** song song dữ liệu (data parallelism).
- **Phân rã:** domain/data decomposition — **1D slab** theo trục x; rank r giữ lát [r·Lx/P,(r+1)·Lx/P).
- **Giao tiếp:** `MPI_Sendrecv` **blocking**, topology **ring/tuyến tính** (2 hàng xóm trái/phải);
  ghost exchange **2 pha** (vị trí cho mật độ, rồi ρ,p cho lực); output **master-gather** (`MPI_Gatherv` về rank 0).
- **Cân bằng tải:** 1D-slab cân bằng khi hạt phân bố đều theo x → dùng kịch bản pool.
- **Lưu ý vật lý:** gia tốc = F/ρ (Müller 2003, đại lượng là mật độ lực), KHÔNG phải F/m.
