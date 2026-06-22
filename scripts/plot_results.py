#!/usr/bin/env python3
"""plot_results.py — Vẽ 3 biểu đồ báo cáo từ CSV do run_experiments.sh sinh ra.

    python3 scripts/plot_results.py            # vẽ tất cả biểu đồ có dữ liệu
    python3 scripts/plot_results.py tn1        # chỉ vẽ 1 biểu đồ

Sinh ra (trong results/):
    tn1_runtime_vs_N.png    — thời gian chạy theo N (có/không truyền thông)
    tn2_load_balance.png    — stacked bar compute+comm mỗi rank
    tn3_speedup.png         — runtime vs P  +  speedup vs P (kèm đường lý tưởng)
"""
import sys, os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "results"

def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

# ---------- TN1: runtime vs N ----------
def plot_tn1():
    p = os.path.join(R, "tn1_runtime_vs_N.csv")
    if not os.path.exists(p): return print("[tn1] thieu", p)
    rows = read_csv(p)
    N   = [int(float(r["N"]))      for r in rows]
    Tt  = [float(r["T_total"])     for r in rows]
    Tc  = [float(r["T_compute"])   for r in rows]
    plt.figure(figsize=(8, 5))
    plt.plot(N, Tt, "o-", label="Co truyen thong (T_total)")
    plt.plot(N, Tc, "s--", label="Khong truyen thong (T_compute)")
    plt.axhspan(120, 180, color="green", alpha=0.12, label="Vung muc tieu 2-3 phut")
    plt.xlabel("N (so hat)"); plt.ylabel("Thoi gian chay (s)")
    plt.title("TN1 — Thoi gian chay theo kich thuoc N"); plt.legend(); plt.grid(alpha=.3)
    out = os.path.join(R, "tn1_runtime_vs_N.png"); plt.tight_layout(); plt.savefig(out, dpi=120)
    print("[tn1] da luu", out)

# ---------- TN2: load balance stacked bar ----------
def plot_tn2():
    p = os.path.join(R, "tn2_timing.csv")
    if not os.path.exists(p):
        p = os.path.join("output", "timing.csv")   # dự phòng
    if not os.path.exists(p): return print("[tn2] thieu tn2_timing.csv")
    rows = read_csv(p)
    ranks = [int(r["rank"])       for r in rows]
    comp  = [float(r["compute_s"]) for r in rows]
    comm  = [float(r["comm_s"])    for r in rows]
    plt.figure(figsize=(max(7, len(ranks)*0.7), 5))
    plt.bar(ranks, comp, label="Compute", color="#2e7d32")
    plt.bar(ranks, comm, bottom=comp, label="Communication", color="#c62828")
    plt.xlabel("Rank (tien trinh)"); plt.ylabel("Thoi gian (s)")
    cmax, cmin = max(comp), min(comp)
    imb = 100*(cmax-cmin)/cmax if cmax else 0
    plt.title(f"TN2 — Can bang tai (mat can bang compute = {imb:.1f}%)")
    plt.xticks(ranks); plt.legend(); plt.grid(axis="y", alpha=.3)
    out = os.path.join(R, "tn2_load_balance.png"); plt.tight_layout(); plt.savefig(out, dpi=120)
    print("[tn2] da luu", out, f"(mat can bang {imb:.1f}%)")

# ---------- TN3: runtime + speedup vs P ----------
def plot_tn3():
    p = os.path.join(R, "tn3_speedup.csv")
    if not os.path.exists(p): return print("[tn3] thieu", p)
    rows = read_csv(p)
    P  = [int(float(r["P"]))      for r in rows]
    Tt = [float(r["T_total"])     for r in rows]
    Tc = [float(r["T_compute"])   for r in rows]
    t1 = Tt[0]                              # T(1) lam moc speedup
    sp = [t1/t for t in Tt]
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].plot(P, Tt, "o-", label="Co truyen thong")
    ax[0].plot(P, Tc, "s--", label="Khong truyen thong")
    ax[0].set_xlabel("So tien trinh P"); ax[0].set_ylabel("Thoi gian (s)")
    ax[0].set_title("Runtime vs P"); ax[0].set_xscale("log", base=2)
    ax[0].set_xticks(P); ax[0].get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].plot(P, sp, "o-", label="Speedup thuc te")
    ax[1].plot(P, P,  "k:", label="Ly tuong (tuyen tinh)")
    ax[1].set_xlabel("So tien trinh P"); ax[1].set_ylabel("Speedup = T(1)/T(P)")
    ax[1].set_title("Speedup vs P"); ax[1].set_xscale("log", base=2)
    ax[1].set_xticks(P); ax[1].get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax[1].legend(); ax[1].grid(alpha=.3)
    out = os.path.join(R, "tn3_speedup.png"); plt.tight_layout(); plt.savefig(out, dpi=120)
    print("[tn3] da luu", out)

if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "tn1"): plot_tn1()
    if which in ("all", "tn2"): plot_tn2()
    if which in ("all", "tn3"): plot_tn3()
