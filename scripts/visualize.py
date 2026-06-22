#!/usr/bin/env python3
"""visualize.py — Dựng hoạt ảnh 3D từ các file output/frame_*.csv.

Dùng:
    python3 scripts/visualize.py            # lưu animation .gif
    python3 scripts/visualize.py --frame 10 # chỉ xem 1 frame -> PNG
    python3 scripts/visualize.py --color rank   # tô màu theo rank/máy (mặc định: vận tốc)

Chỉ cần numpy + matplotlib (không cần ParaView/pyvista).
Với cluster lớn, dùng ParaView đọc trực tiếp các CSV để render đẹp hơn.
"""
import argparse, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

OUT = "output"

def load(path):
    # x,y,z,vx,vy,vz,rho,p,rank
    return np.loadtxt(path, delimiter=",", skiprows=1)

def color_values(data, mode):
    if mode == "rank":
        return data[:, 8]
    speed = np.linalg.norm(data[:, 3:6], axis=1)  # |v|
    return speed

def setup_ax(ax, L=(1.0, 0.6, 0.6)):
    ax.set_xlim(0, L[0]); ax.set_ylim(0, L[1]); ax.set_zlim(0, L[2])
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    try: ax.set_box_aspect(L)
    except Exception: pass

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", type=int, default=None)
    ap.add_argument("--color", choices=["speed", "rank"], default="speed")
    ap.add_argument("--out", default="sph3d.gif")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(OUT, "frame_*.csv")))
    if not files:
        print("Khong tim thay frame nao trong output/. Hay chay mo phong truoc.")
        return
    print(f"Tim thay {len(files)} frame.")

    if args.frame is not None:
        sel = args.frame
        files = [files[sel]]
    else:
        sel = 0

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    def draw(i):
        ax.clear(); setup_ax(ax)
        d = load(files[i])
        if d.ndim == 1: d = d.reshape(1, -1)
        c = color_values(d, args.color)
        ax.scatter(d[:, 0], d[:, 1], d[:, 2], c=c, cmap="viridis", s=6)
        label = sel + i if len(files) == 1 else i
        ax.set_title(f"SPH 3D Dam Break — frame {label} ({len(d)} hat)")

    if len(files) == 1:
        draw(0)
        png = "frame.png"; fig.savefig(png, dpi=110)
        print(f"Da luu {png}")
    else:
        anim = FuncAnimation(fig, draw, frames=len(files), interval=80)
        anim.save(args.out, writer=PillowWriter(fps=12))
        print(f"Da luu {args.out}")

if __name__ == "__main__":
    main()
