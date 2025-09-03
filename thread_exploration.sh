#!/bin/bash
# === Load SLURM options ===
source my_slurm.sh 

make clean
make -j$(nproc) OPENMP_SCHEDULE=3

# Submit the thread scaling study
echo "Submitting thread scalability study..."
sbatch \
  $COMMON_OPTS \
  $EXTRA_OPTS \
  --job-name=thread_study \
  --nodes=1 \
  --cpus-per-task=112 \
  --ntasks-per-node=1 \
  --mem=0 \
  --export=ENV=$ENV \
  --output=thread_scaling_%j.out \
  thread_study.sh

echo "Thread study submitted. Results will be in thread_results.csv and thread_speedup.csv"
