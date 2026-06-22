#!/usr/bin/env bash
# run_experiments.sh — Tự động chạy các thí nghiệm TN0–TN3 và ghi CSV cho plot_results.py
#
# Cấu hình qua biến môi trường:
#   CORES     : tổng số nhân vật lý của cụm (vd 3 máy x 4 nhân = 12). Mặc định: nproc máy hiện tại.
#   HOSTFILE  : đường dẫn hostfile MPI (vd ~/hostfile). Để trống nếu chạy 1 máy.
#   STEPS     : số bước thời gian. Mặc định 1000.
#
# Cách dùng:
#   ./scripts/run_experiments.sh tn0 [dx] [dt]        # kiểm tra đúng đắn serial vs MPI
#   ./scripts/run_experiments.sh tn1                  # quét N -> tìm N cho 2-3 phút
#   ./scripts/run_experiments.sh tn2 <dx> <dt>        # cân bằng tải tại N đã chọn (P=CORES)
#   ./scripts/run_experiments.sh tn3 <dx_2N> <dt>     # speedup: input 2N, quét P
set -e
cd "$(dirname "$0")/.."          # về thư mục dự án
mkdir -p results output

CORES=${CORES:-$(nproc)}
STEPS=${STEPS:-1000}
HFARG=""; [ -n "$HOSTFILE" ] && HFARG="--hostfile $HOSTFILE"
ROOTARG=""; [ "$(id -u)" = "0" ] && ROOTARG="--allow-run-as-root"
MPIRUN="mpirun $HFARG --oversubscribe $ROOTARG"

# In tham số mpi: $1=np $2=dx $3=dt $4=scenario -> trả về dòng kết quả [mpi]
run_mpi() { $MPIRUN -np "$1" ./sph_mpi "$STEPS" "$2" "$3" "$4"; }

# Lấy số sau dấu ':' ở dòng chứa từ khóa
field() { grep "$1" | sed 's/.*: *//; s/ s.*//' | head -1; }

ensure_build() { [ -x ./sph_mpi ] && [ -x ./sph_serial ] || make >/dev/null; }

case "${1:-}" in

  tn0)  # --- KIỂM TRA ĐÚNG ĐẮN ---
    ensure_build
    DX=${2:-0.02}; DT=${3:-0.0002}
    echo "## TN0 Correctness | dx=$DX dt=$DT steps=$STEPS scenario=pool"
    echo "Phien ban,P,cx,cy,cz"
    rm -f output/frame_*.csv
    ./sph_serial "$STEPS" "$DX" "$DT" pool >/dev/null 2>&1
    F=$(ls -1 output/frame_*.csv | tail -1)
    awk -F, 'NR>1{x+=$1;y+=$2;z+=$3;n++}END{printf "serial,1,%.4f,%.4f,%.4f\n",x/n,y/n,z/n}' "$F"
    for P in 2 $CORES; do
      rm -f output/frame_*.csv
      run_mpi "$P" "$DX" "$DT" pool >/dev/null 2>&1
      F=$(ls -1 output/frame_*.csv | tail -1)
      awk -F, -v p="$P" 'NR>1{x+=$1;y+=$2;z+=$3;n++}END{printf "mpi,%d,%.4f,%.4f,%.4f\n",p,x/n,y/n,z/n}' "$F"
    done
    echo ">> Trong tam phai trung nhau (~<0.5%) => song song dung."
    ;;

  tn1)  # --- QUÉT N, runtime vs N, tại P=CORES ---
    ensure_build
    echo "## TN1 runtime vs N | P=$CORES steps=$STEPS scenario=pool"
    echo "N,T_total,T_comm,T_compute" > results/tn1_runtime_vs_N.csv
    # cặp (dx, dt=0.01*dx). Bỏ bớt dòng cuối nếu máy yếu.
    for pair in "0.025 0.00025" "0.020 0.00020" "0.018 0.00018" "0.015 0.00015" \
                "0.012 0.00012" "0.010 0.00010" "0.008 0.00008"; do
      set -- $pair; DX=$1; DT=$2
      OUT=$(run_mpi "$CORES" "$DX" "$DT" pool)
      N=$(echo "$OUT"   | grep "N_tong" | grep -oP 'N_tong=\K[0-9]+' | head -1)
      TT=$(echo "$OUT"  | field "Tong thoi gian")
      TM=$(echo "$OUT"  | field "Communication (max")
      TC=$(echo "$OUT"  | field "Compute (khong")
      echo "$N,$TT,$TM,$TC" | tee -a results/tn1_runtime_vs_N.csv
    done
    echo ">> Chon N (dx) cho T_total ~ 120-180s lam N* cho TN2/TN3."
    ;;

  tn2)  # --- CÂN BẰNG TẢI tại N*, P=CORES ---
    ensure_build
    DX=${2:?Can dx}; DT=${3:?Can dt}
    echo "## TN2 load balance | N(dx=$DX) P=$CORES steps=$STEPS scenario=pool"
    run_mpi "$CORES" "$DX" "$DT" pool | grep -E "N_tong|per-rank|Mat can bang|^ +[0-9]"
    cp output/timing.csv results/tn2_timing.csv
    echo ">> Bang per-rank da ghi vao results/tn2_timing.csv (de ve stacked bar)."
    ;;

  tn3)  # --- SPEEDUP: input 2N, quét P=1,2,4,8,...,2*CORES ---
    ensure_build
    DX=${2:?Can dx cho 2N}; DT=${3:?Can dt}
    echo "## TN3 speedup | 2N(dx=$DX) steps=$STEPS scenario=pool"
    echo "P,T_total,T_compute" > results/tn3_speedup.csv
    P=1
    while [ "$P" -le $((2*CORES)) ]; do
      OUT=$(run_mpi "$P" "$DX" "$DT" pool)
      TT=$(echo "$OUT" | field "Tong thoi gian")
      TC=$(echo "$OUT" | field "Compute (khong")
      echo "$P,$TT,$TC" | tee -a results/tn3_speedup.csv
      P=$((P*2))
    done
    echo ">> Da ghi results/tn3_speedup.csv. Chay plot_results.py de ve."
    ;;

  *) echo "Dung: $0 [tn0|tn1|tn2 <dx> <dt>|tn3 <dx> <dt>]"; exit 1;;
esac
