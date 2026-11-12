#!/bin/bash
#SBATCH --job-name=pw_scan2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=120:00:00
#SBATCH --output=slurm-%A_%a.out
# Scan:      sbatch --array=1-34%5 slurm_array.sh          (priority 1: in-plane lines + saddle region)
#            sbatch --array=35-143%5 slurm_array.sh        (the rest)
# Restarts:  sbatch --array=1-N%5 slurm_array.sh restart   (after python3 make_restart.py)
# Step 2b (tau frozen): sbatch --array=1-65%5 slurm_array.sh tau   then  --array=66-142%5 ... tau
#            restarts:  sbatch --array=1-N%5 slurm_array.sh taurestart
# Step 3 (3D grid, cluster share of 80 jobs): sbatch --array=1-80%2 slurm_array.sh 3d
# No SLURM:  bash slurm_array.sh scan 7      or   bash slurm_array.sh restart 2

# Gaussian: works for both accounts (ramachandran.g: ~/test/g09 ; balamurali.mm: ~/Programs/g09)
if [ -d $HOME/test/g09 ]; then export g09root=$HOME/test/; else export g09root=$HOME/Programs/; fi
export GAUSS_EXEDIR=$g09root/g09
export PATH=$GAUSS_EXEDIR:$PATH
GAUSS=g09

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
mode=${1:-scan}
case "$mode" in
  restart)    list=joblist_restart.txt;            dir=inputs_restart ;;
  tau)        list=joblist_tau.txt;                dir=inputs_tau ;;
  taurestart) list=joblist_inputs_tau_restart.txt; dir=inputs_tau_restart ;;
  3d)         list=joblist_3d_cluster.txt;         dir=inputs_3d ;;
  3drestart)  list=joblist_inputs_3d_restart.txt;  dir=inputs_3d_restart ;;
  *)          list=joblist.txt;                    dir=inputs ;;
esac
id=${SLURM_ARRAY_TASK_ID:-$2}
[ -z "$id" ] && { echo "No job number. Use: sbatch --array=1-34%5 slurm_array.sh   or   bash slurm_array.sh scan N"; exit 1; }
name=$(sed -n "${id}p" $list | awk '{print $1}')
[ -z "$name" ] && { echo "No job $id in $list"; exit 1; }

export GAUSS_SCRDIR="/tmp/g09_${USER}_${SLURM_JOB_ID:-local}_${mode}_${id}"
mkdir -p "$GAUSS_SCRDIR"
cd $dir || exit 1
if [ -f "$name.log" ] && grep -q "Normal termination" "$name.log" && grep -q "Optimization completed" "$name.log"; then echo "$name already done"; exit 0; fi
echo "$(date '+%F %T')  running $name on $(hostname)"
$GAUSS < "$name.gjf" > "$name.log"
echo "$(date '+%F %T')  finished $name"
case "$GAUSS_SCRDIR" in */g09_*) rm -rf "$GAUSS_SCRDIR";; esac
