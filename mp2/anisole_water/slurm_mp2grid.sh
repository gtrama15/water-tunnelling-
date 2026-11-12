#!/bin/bash
#SBATCH --job-name=mp2anis
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=120:00:00
#SBATCH --output=slurm-mp2anis-%A_%a.out
# sbatch --array=1-N%3 slurm_mp2corr.sh       (names from joblist.txt)
if [ -d $HOME/test/g09 ]; then export g09root=$HOME/test/; else export g09root=$HOME/Programs/; fi
export GAUSS_EXEDIR=$g09root/g09; export PATH=$GAUSS_EXEDIR:$PATH
cd "$SLURM_SUBMIT_DIR" || exit 1
n=$(sed -n "${SLURM_ARRAY_TASK_ID}p" joblist.txt)
[ -z "$n" ] && exit 1
cd inputs || exit 1
if [ -f $n.log ] && grep -q "Normal termination" $n.log; then echo "$n already done"; exit 0; fi
export GAUSS_SCRDIR=$SLURM_SUBMIT_DIR/scr_$n; mkdir -p "$GAUSS_SCRDIR"     # node /tmp is too small for MP2
echo "$(date '+%F %T')  $n on $(hostname)"; df -h "$GAUSS_SCRDIR" | tail -1
g09 < $n.gjf > $n.log
echo "$(date '+%F %T')  done $n  $(grep -c 'Normal termination' $n.log)"
rm -rf "$GAUSS_SCRDIR"
