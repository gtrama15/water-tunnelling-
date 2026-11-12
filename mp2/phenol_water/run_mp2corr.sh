#!/bin/bash
# Runs the MP2 correction points one after another (each on 12 cores, 32 GB, scratch in /tmp).
# Usage:  nohup bash run_mp2corr.sh > run_mp2corr.log 2>&1 &      (skips finished jobs; safe to restart)
export g09root=${g09root:-/home/rama/Programs/G09-D/}
. $g09root/g09/bsd/g09.profile
cd "$(dirname "$0")" || exit 1
for n in $(cat joblist_mp2corr.txt); do
  cd inputs
  if [ -f $n.log ] && grep -q "Normal termination" $n.log; then cd ..; continue; fi
  export GAUSS_SCRDIR=${SCR:-/tmp}/g09_${USER}_$n; mkdir -p $GAUSS_SCRDIR
  echo "$(date '+%F %T')  start $n"; g09 < $n.gjf > $n.log; echo "$(date '+%F %T')  done  $n  $(grep -c 'Normal termination' $n.log)"
  rm -rf $GAUSS_SCRDIR; cd ..
done
