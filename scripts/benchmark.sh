#!/usr/bin/env bash
# benchmark.sh — Đo hiệu năng multi-node (mục 9). Chạy từ thư mục dự án.
#
#   ./scripts/benchmark.sh strong   # strong scaling: cố định N, tăng P
#   ./scripts/benchmark.sh weak     # weak scaling: cố định N/rank
#
# Trên cluster, thêm --hostfile ~/hostfile vào mỗi lệnh mpirun bên dưới.
set -e
HF=${HOSTFILE:-}          # export HOSTFILE=~/hostfile để chạy trên cụm
HFARG=""; [ -n "$HF" ] && HFARG="--hostfile $HF"
STEPS=${STEPS:-300}

run() {  # run <np> <dx> <nhãn>
  echo "=== $3 : np=$1 dx=$2 ==="
  mpirun $HFARG --oversubscribe -np "$1" ./sph_mpi "$STEPS" "$2" \
    | grep -E "Tong thoi gian|Compute|Communication|Comm/Total|N_tong"
  echo
}

case "${1:-strong}" in
  strong)
    echo "## STRONG SCALING — N co dinh (~15K hat, dx=0.02), tang P"
    for np in 1 2 3 6; do run "$np" 0.02 "strong P=$np"; done
    ;;
  weak)
    echo "## WEAK SCALING — ~N/rank co dinh"
    # dx nho hon -> nhieu hat hon, tang theo P de N/rank ~ hang so
    run 1 0.030 "weak P=1"
    run 3 0.021 "weak P=3"
    run 6 0.017 "weak P=6"
    ;;
  *) echo "Dung: $0 [strong|weak]"; exit 1;;
esac
