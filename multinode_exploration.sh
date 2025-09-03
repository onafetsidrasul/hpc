#!/bin/bash
# === Load common SLURM options ===
source my_slurm.sh  # defines COMMON_OPTS
make clean 
make -j$(nproc) OPENMP_SCHEDULE=1 

#
# === Read best_n_threads and best_n_tasks ===
if [[ ! -f best_thread_task ]]; then
  echo "❌ Missing 'best_thread_task' file" >&2
  exit 1
fi

read BEST_THREADS BEST_TASKS < <(awk -F, 'NR==1{next} {print $1, $2}' best_thread_task)
echo "✅ Read best config: $BEST_TASKS MPI tasks × $BEST_THREADS threads"

# === Submit job array for multinode scaling ===
NODE_COUNTS=(1 2 4 8 16)
# NODE_COUNTS=(2 8 16)
for i in "${!NODE_COUNTS[@]}"; do
  NODES="${NODE_COUNTS[$i]}"
  # NTASKS=$((NODES*BEST_TASKS))
  echo "Starting job with $NODES nodes and $BEST_TASKS tasks per node, with $BEST_THREADS cpus per task..."
  sbatch \
    $COMMON_OPTS \
    $EXTRA_OPTS \
    --job-name=multinode_study\
    --nodes=$NODES \
    --ntasks-per-node=$BEST_TASKS \
    --cpus-per-task=$BEST_THREADS \
    --export=ENV=$ENV,BEST_THREADS=$BEST_THREADS,BEST_TASKS=$BEST_TASKS,NODES=$NODES \
    --output=log_multinode_${NODES}nodes.out \
    multinode_study.sh
done

