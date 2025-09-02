#!/bin/bash
set -x

# Create a timestamp marker
NOW=$(date +%Y-%m-%d_%H%M%S)

DIR=$(pwd)/multinode_exp
if [ ! -d "$DIR" ]; then
  mkdir -p "$DIR" 
fi 

# Initialize environment and log file
LOG_FILE="$DIR/out_${NODES}nodes${NOW}.txt"
ERROR_FILE="$DIR/err_${NODES}nodes${NOW}.txt"
STRONG_TIME_FILE="$DIR/time_${NODES}_STRONG${NOW}.txt"
WEAK_TIME_FILE="$DIR/time_${NODES}_WEAK${NOW}.txt"

STRONG_RESULTS="$DIR/multinode_results_STRONG${NOW}.csv"
WEAK_RESULTS="$DIR/multinode_results_WEAK${NOW}.csv"

## Label the files that will later be used to compute speedup and efficiency
if [[ $NODES -eq 1 ]]; then
	echo "nodes,x,y, time, elapsed_time" > "$STRONG_RESULTS" 
	echo "nodes,x,y, time, elapsed_time" > "$WEAK_RESULTS" 
else
    # Ensure no stale data exists for this NODES value
    grep -v "^$NODES," "$STRONG_RESULTS" > "${STRONG_RESULTS}.tmp" && mv "${STRONG_RESULTS}.tmp" "$STRONG_RESULTS"
    grep -v "^$NODES," "$WEAK_RESULTS" > "${WEAK_RESULTS}.tmp" && mv "${WEAK_RESULTS}.tmp" "$WEAK_RESULTS"
fi

# Clear previous log files
> "$LOG_FILE"
> "$ERROR_FILE"

# Set environment
if [[ $ENV == "LEONARDO" ]]; then
    echo "[$(date)] Setting environment for LEONARDO" >> "$LOG_FILE"
    module load openmpi/4.1.6--gcc--12.2.0
elif [[ $ENV == "ORFEO" ]]; then
    echo "[$(date)] Setting environment for ORFEO" >> "$LOG_FILE"
    module load openMPI
fi

EXEC=./stencil_parallel

# Common parameters
X_BASE=10000
Y_BASE=10000
total_tasks=$((NODES * BEST_TASKS))

# Set OpenMP parameters
export OMP_NUM_THREADS=$BEST_THREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "=============================================================" >> "$LOG_FILE"
echo "[$(date)] Starting experiment with $NODES nodes" >> "$LOG_FILE"
echo "Configuration:" >> "$LOG_FILE"
echo "  - MPI tasks per node: $BEST_TASKS" >> "$LOG_FILE"
echo "  - Threads per task: $BEST_THREADS" >> "$LOG_FILE"
echo "  - Total MPI tasks: $total_tasks" >> "$LOG_FILE"
echo "=============================================================" >> "$LOG_FILE"

# Strong scaling test (fixed problem size, increasing resources)
echo "[$(date)] Starting STRONG SCALING test" >> "$LOG_FILE"
x=$X_BASE
y=$Y_BASE

echo "  Problem size: ${x}x${y} (fixed)" >> "$LOG_FILE"
echo "  Running command: mpirun -np $total_tasks $EXEC -x $x -y $y -o 1" >> "$LOG_FILE"

start_time=$(date +%s.%N)
/usr/bin/time -f "%e" -o "$STRONG_TIME_FILE" mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 >> "$LOG_FILE" 2>> "$ERROR_FILE"
strong_exit_code=$?
end_time=$(date +%s.%N)
elapsed=$(echo "$end_time - $start_time" | bc)

if [[ $strong_exit_code -eq 0 ]]; then
    strong_time=$(cat "$STRONG_TIME_FILE")
    echo "[$(date)] STRONG SCALING completed successfully in ${strong_time}s (${elapsed}s wall time)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,$strong_time,$elapsed" >> "$STRONG_RESULTS"
else
    echo "[$(date)] ERROR in STRONG SCALING test (exit code: $strong_exit_code)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,ERROR,$strong_exit_code" >> "$STRONG_RESULTS"
fi

x=$(( (X_BASE)* $NODES)) 
y=$Y_BASE


# Weak scaling test (problem size grows with resources)
echo "[$(date)] Starting WEAK SCALING test" >> "$LOG_FILE"

echo "  Problem size: ${x}x${y} (scales with nodes)" >> "$LOG_FILE"
echo "  Running command: mpirun -np $total_tasks $EXEC -x $x -y $y -o 1" >> "$LOG_FILE"

start_time=$(date +%s.%N)
/usr/bin/time -f "%e" -o "$WEAK_TIME_FILE" mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 >> "$LOG_FILE" 2>> "$ERROR_FILE"
weak_exit_code=$?
end_time=$(date +%s.%N)
elapsed=$(echo "$end_time - $start_time" | bc)

if [[ $weak_exit_code -eq 0 ]]; then
    weak_time=$(cat "$WEAK_TIME_FILE")
    echo "[$(date)] WEAK SCALING completed successfully in ${weak_time}s (${elapsed}s wall time)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,$weak_time,$elapsed" >> "$WEAK_RESULTS"
else
    echo "[$(date)] ERROR in WEAK SCALING test (exit code: $weak_exit_code)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,ERROR,$weak_exit_code" >> "$WEAK_RESULTS"
fi

# Final summary
echo "=============================================================" >> "$LOG_FILE"
echo "[$(date)] Experiment completed for $NODES nodes" >> "$LOG_FILE"
echo "Results summary:" >> "$LOG_FILE"
echo "  - Strong scaling: ${strong_time:-ERROR}s" >> "$LOG_FILE"
echo "  - Weak scaling: ${weak_time:-ERROR}s" >> "$LOG_FILE"
./process_results
echo "Processed results(i.e speedup and efficiency) can be found in the _metrics file in this directory."
echo "=============================================================" >> "$LOG_FILE"

