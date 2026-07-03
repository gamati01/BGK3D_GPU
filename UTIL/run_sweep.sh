#!/usr/bin/env bash
# ====================================================================
#  run_sweep.sh
#
#  Drives a size sweep of the BGK3D solver, running each configuration
#  NRUNS times and reporting the mean of:
#     - Mlups                       (global, from the "# Mlups" line)
#     - GPU Mean Power (W)          (of the busiest GPU)
#     - energy per gridpoint (nJ)   (of the busiest GPU)
#
#  Mlups is read from the per-run log (the solver always prints the
#  "# Mlups" line to stdout).  Power and energy/gridpoint are read from
#  bgk.perf.  Because a node can host several GPUs, for every run we pick
#  the GPU that consumed the MOST energy (J) and report that GPU's power
#  and energy/gridpoint.  (Build the solver with ENERGY=1 so the energy
#  lines are present.)
#
#  bgk.input (a Fortran namelist) is rewritten for each configuration:
#  Size sets lx = ly = lz; itfin, icheck and svisc vary per the table.
#
#  Usage:
#     cd RUN
#     ../UTIL/run_sweep.sh [executable] [nruns]
#
#  Examples:
#     ../UTIL/run_sweep.sh ./bgk3d.hip.x        # 5 runs (default)
#     ../UTIL/run_sweep.sh ./bgk3d.cuda.x 3     # 3 runs
#
#  Environment overrides:
#     RUNDIR   directory holding bgk.input / the executable (default: $PWD)
#     NRUNS    number of repetitions per configuration       (default: 5)
#     OUTCSV   summary CSV file name                          (default: sweep_results.csv)
# ====================================================================
set -euo pipefail

# -------- knobs ------------------------------------------------------
EXE="${1:-${EXE:-./bgk3d.hip.x}}"
NRUNS="${2:-${NRUNS:-5}}"
RUNDIR="${RUNDIR:-$PWD}"
OUTCSV="${OUTCSV:-sweep_results_hip.csv}"

INPUT="bgk.input"
PERF="bgk.perf"

# namelist parameters kept constant across the whole sweep
U0=0.1
IVTIM=500000
ISIGNAL=250
ITSAVE=500000
IRESTART=0
INIT_V=0
IPAD=0

# configurations:  "size itfin icheck svisc"
configs=(
  "128  200000 40000 0.0128"
  "180  200000 40000 0.0180"
  "256  100000 20000 0.0256"
  "360  100000 20000 0.0360"
  "512  50000  10000 0.0512"
  "720  50000  10000 0.0720"
  "1000 25000  10000 0.1000"
)

# -------- setup ------------------------------------------------------
cd "$RUNDIR"

if [[ ! -x "$EXE" && ! -f "$EXE" ]]; then
  echo "ERROR: executable '$EXE' not found in $RUNDIR" >&2
  exit 1
fi

LOGDIR="${LOGDIR:-sweep_logs_hip}"
mkdir -p "$LOGDIR"

# preserve any existing bgk.input and restore it on exit
if [[ -f "$INPUT" ]]; then
  cp -f "$INPUT" "$LOGDIR/bgk.input.orig"
  trap 'cp -f "$LOGDIR/bgk.input.orig" "$INPUT" 2>/dev/null || true' EXIT
fi

echo "size,nruns_ok,mlups_mean,power_mean_W,energy_per_gp_mean_nJ,busiest_gpu" > "$OUTCSV"

# write the namelist for a given configuration
write_input() {
  local size=$1 itfin=$2 icheck=$3 svisc=$4
  cat > "$INPUT" <<EOF
&parameters
lx = $size
ly = $size
lz = $size
svisc=$svisc
u0=$U0
itfin=$itfin
ivtim=$IVTIM
isignal=$ISIGNAL
itsave=$ITSAVE
icheck=$icheck
irestart=$IRESTART
init_v=$INIT_V
ipad=$IPAD        /
EOF
}

# parse run log + bgk.perf -> "mlups power_busiest epg_busiest gpu_idx"
# Mlups is read from the run log because the solver always prints the
# "# Mlups" line to stdout, whereas bgk.perf only carries it when built
# with -DPROFILING (ENERGY-only builds omit it).  Power and energy come
# from bgk.perf.  busiest = GPU with the largest energy(J).
# Prints "NA" fields if missing.
parse_perf() {
  local perf=$1 log=$2
  awk '
    FNR==NR { if ($0 ~ /Mlups/) mlups = $NF; next }   # first file: the run log
    /Mean Power/           { g=$2; gsub(/[^0-9]/,"",g); en[g]=$(NF-1)+0; pw[g]=$NF+0 }
    /per gridpoint/        { g=$2; gsub(/[^0-9]/,"",g); epg[g]=$NF+0 }
    END {
      best=""; beste=-1
      for (g in en) if (en[g] > beste) { beste=en[g]; best=g }
      m   = (mlups=="" ? "NA" : mlups)
      if (best=="") { print m, "NA", "NA", "NA"; exit }
      print m, pw[best], epg[best], best
    }
  ' "$log" "$perf"
}

# -------- main sweep -------------------------------------------------
for cfg in "${configs[@]}"; do
  read -r size itfin icheck svisc <<< "$cfg"
  echo "=================================================================" >&2
  echo ">> size=$size  itfin=$itfin  icheck=$icheck  svisc=$svisc  (x$NRUNS)" >&2

  write_input "$size" "$itfin" "$icheck" "$svisc"

  tmp="$(mktemp)"
  for ((r=1; r<=NRUNS; r++)); do
    rm -f "$PERF"
    log="$LOGDIR/run_${size}_${r}.log"
    echo "   run $r/$NRUNS ..." >&2
    if ! "$EXE" > "$log" 2>&1; then
      echo "   !! run $r failed (see $log)" >&2
      continue
    fi
    if [[ ! -f "$PERF" ]]; then
      echo "   !! no $PERF produced on run $r (see $log)" >&2
      continue
    fi
    cp -f "$PERF" "$LOGDIR/bgk.perf.${size}.${r}"
    vals="$(parse_perf "$PERF" "$log")"
    echo "      -> mlups,power,epg,gpu = $vals" >&2
    echo "$vals" >> "$tmp"
  done

  # mean over the successful runs of this configuration
  awk -v sz="$size" '
    $1!="NA" { m+=$1; nm++ }
    $2!="NA" { p+=$2; np++ }
    $3!="NA" { e+=$3; ne++ }
    $4!="NA" { gpu=$4 }
    END {
      mm  = (nm? sprintf("%.3f", m/nm) : "NA")
      pp  = (np? sprintf("%.3f", p/np) : "NA")
      ee  = (ne? sprintf("%.3f", e/ne) : "NA")
      g   = (gpu=="" ? "NA" : gpu)
      printf "%s,%d,%s,%s,%s,%s\n", sz, nm, mm, pp, ee, g
    }
  ' "$tmp" | tee -a "$OUTCSV" >&2
  rm -f "$tmp"
done

echo "=================================================================" >&2
echo "Summary written to $RUNDIR/$OUTCSV" >&2
echo "Per-run logs and bgk.perf copies are in $RUNDIR/$LOGDIR/" >&2
column -s, -t "$OUTCSV" >&2 || cat "$OUTCSV" >&2
