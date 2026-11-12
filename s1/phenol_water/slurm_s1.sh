#!/bin/bash
#SBATCH --job-name=s1vert
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-s1-%A_%a.out
# test first:  sbatch --array=1 slurm_s1.sh        then:  sbatch --array=2-171%2 slurm_s1.sh   (2 x 24 cores)
# if the QOS allows 6 running jobs:  sbatch -c 8 --mem=12G --array=2-171%6 slurm_s1.sh       (6 x 8 cores)
# Gaussian %nprocshared and %mem are set from the SLURM request (%mem = 75% of --mem)
if [ -d $HOME/test/g09 ]; then export g09root=$HOME/test/; else export g09root=$HOME/Programs/; fi
export GAUSS_EXEDIR=$g09root/g09; export PATH=$GAUSS_EXEDIR:$PATH
cd "$SLURM_SUBMIT_DIR" || exit 1
n=$(sed -n "${SLURM_ARRAY_TASK_ID}p" joblist_cluster.txt); [ -z "$n" ] && exit 1
cd inputs || exit 1
if [ -f $n.log ] && grep -q "Normal termination" $n.log; then echo "$n already done"; exit 0; fi
export GAUSS_SCRDIR=/tmp/g09_${USER}_${SLURM_JOB_ID}_$n; mkdir -p $GAUSS_SCRDIR
echo "$(date '+%F %T')  $n on $(hostname)"
MEMG=$(( ${SLURM_MEM_PER_NODE:-12288} * 3 / 4 / 1024 ))
sed -e "s/^%nprocshared=.*/%nprocshared=${SLURM_CPUS_PER_TASK:-8}/" -e "s/^%mem=.*/%mem=${MEMG}GB/" $n.gjf | g09 > $n.log
echo "$(date '+%F %T')  done $n  $(grep -c 'Normal termination' $n.log)"
rm -rf $GAUSS_SCRDIR
