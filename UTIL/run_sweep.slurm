#!/bin/bash

#SBATCH --job-name=BGK3D_sweep_hip
#SBATCH --partition=MI300x
#SBATCH --nodes=1
#SBATCH --exclusive
##SBATCH --exclude=TheraC16,TheraC55
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output=BGK3D_sweep_hip.out
#SBATCH --error=BGK3D_sweep_hip.err

source /home/mbycklin/bin/spack-setup-env.sh

module load rocm-afar-drop/23.2.1-gfx942
#module load rocm-afar-drop/23.2.1-gfx950
module load rocm-openmpi openblas rocm-hdf5-ompi-5.0.10/1.14.6

export OMP_NUM_THREADS=1
export HSA_ENABLE_IPC_MODE_LEGACY=1
export HSA_ENABLE_SDMA=1

echo "SLURM Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Current working directory is: $(pwd)"

# ---- sweep configuration -------------------------------------------------
# Directory holding bgk.input and the executable.
RUNDIR="${RUNDIR:-$SLURM_SUBMIT_DIR}"
# Solver executable (built with ENERGY=1 so power/energy lines are present).
EXE="${EXE:-./bgk3d.hip.x}"
# Repetitions per configuration.
NRUNS="${NRUNS:-5}"
# Summary CSV file name.
OUTCSV="${OUTCSV:-sweep_results_hip.csv}"

# Location of run_sweep.sh.  Under SLURM the batch script is copied to a
# spool directory, so ${BASH_SOURCE[0]} cannot be trusted to point at UTIL/.
# Because we are submitted from RUN/ as "sbatch ../UTIL/run_sweep.slurm",
# UTIL/ is a sibling of the submit directory.
SCRIPT_DIR="$(cd "$SLURM_SUBMIT_DIR/../UTIL" && pwd)"
SWEEP="${SWEEP:-$SCRIPT_DIR/run_sweep.sh}"

export RUNDIR NRUNS OUTCSV

echo "=== BGK3D sweep configuration ==="
echo "RUNDIR : ${RUNDIR}"
echo "EXE    : ${EXE}"
echo "NRUNS  : ${NRUNS}"
echo "OUTCSV : ${OUTCSV}"
echo "SWEEP  : ${SWEEP}"
echo "================================="

"${SWEEP}" "${EXE}" "${NRUNS}"

echo "SLURM Job finished at: $(date)"
