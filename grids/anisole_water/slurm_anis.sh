#!/bin/bash
#SBATCH --job-name=anis3d
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=120:00:00
#SBATCH --output=slurm-anis-%A_%a.out
# sbatch -c 4 --mem=12G --array=1-594%6 slurm_anis.sh            (cores/memory from the sbatch line; %mem = 75% of --mem)
# restarts:  sbatch -c 4 --mem=12G --array=1-N%12 slurm_anis.sh restart
if [ -d $HOME/test/g09 ]; then export g09root=$HOME/test/; else export g09root=$HOME/Programs/; fi
export GAUSS_EXEDIR=$g09root/g09; export PATH=$GAUSS_EXEDIR:$PATH
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
case "${1:-scan}" in
  restart) list=joblist_restart.txt; dir=inputs_restart ;;
  *)       list=joblist.txt;         dir=inputs ;;
esac
n=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $list | awk '{print $1}'); [ -z "$n" ] && { echo "no job"; exit 1; }
cd $dir || exit 1
if [ -f $n.log ] && grep -q "Normal termination" $n.log && grep -q "Optimization completed" $n.log; then echo "$n done"; exit 0; fi
export GAUSS_SCRDIR=/tmp/g09_${USER}_${SLURM_JOB_ID}_$n; mkdir -p $GAUSS_SCRDIR
MEMG=$(( ${SLURM_MEM_PER_NODE:-12288} * 3 / 4 / 1024 ))
echo "$(date '+%F %T')  $n on $(hostname)"
sed -e "s/^%nprocshared=.*/%nprocshared=${SLURM_CPUS_PER_TASK:-4}/" -e "s/^%mem=.*/%mem=${MEMG}GB/" $n.gjf > $GAUSS_SCRDIR/in.gjf
g09 < $GAUSS_SCRDIR/in.gjf > $n.log
echo "$(date '+%F %T')  done $n"
rm -rf $GAUSS_SCRDIR
