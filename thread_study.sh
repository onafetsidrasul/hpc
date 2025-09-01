#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=112
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --job-name=thread_scaling

set -x

# Create a timestamp marker
NOW=$(date +%Y-%m-%d_%H%M%S)

# Base output directory
DIR=$(pwd)/thread_exp_${NOW}
mkdir -p "$DIR"

# Define filenames with timestamp
LOG_FILE="$DIR/thread_scaling_${NOW}.log"
RESULTS_FILE="$DIR/thread_results_${NOW}.csv"
SPEEDUP_FILE="$DIR/thread_speedup_${NOW}.csv"

# Set environment
if [[ $ENV == "LEONARDO" ]]; then
    echo "[$(date)] Setting environment for LEONARDO" >> "$LOG_FILE"
    module load openmpi/4.1.6--gcc--12.2.0
elif [[ $ENV == "ORFEO" ]]; then
    echo "[$(date)] Setting environment for ORFEO" >> "$LOG_FILE"
    module load openMPI
fi

echo "Started thread scalability study at $(date)" > "$LOG_FILE"

EXEC=./stencil_parallel

# Initialize results file
echo "threads,time,wall_time" > "$RESULTS_FILE"

# Test different thread counts
for threads in 1 2 4 8 16 32 56 84 112; do
    echo "[$(date)] Testing OMP_NUM_THREADS=$threads" >> "$LOG_FILE"
    
    export OMP_NUM_THREADS=$threads
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close

    # Run with timing
    start_time=$(date +%s.%N)
    /usr/bin/time -f "%e" -o "$DIR/time_tmp.txt" $EXEC -o 1 > "$DIR/out_${threads}threads_${NOW}.txt" 2> "$DIR/err_${threads}threads_${NOW}.txt"
    runtime=$(cat "$DIR/time_tmp.txt")
    end_time=$(date +%s.%N)
    wall_time=$(echo "$end_time - $start_time" | bc)
    
    # Log results
    echo "$threads,$runtime,$wall_time" >> "$RESULTS_FILE"
    echo "Completed $threads threads in $runtime s (wall: $wall_time s)" >> "$LOG_FILE"
done

# Compute speedup and efficiency
awk -F, -v out="$SPEEDUP_FILE" '
  BEGIN {
    print "threads,time,speedup,efficiency,wall_time" > out
  }
  NR==1 {next}  # Skip header
  NR==2 {
    t1 = $2  # Baseline time (1 thread)
  }
  {
    speedup = t1/$2
    efficiency = speedup/$1
    printf("%d,%.2f,%.2f,%.2f,%.2f\n", $1, $2, speedup, efficiency, $3) >> out
  }
' "$RESULTS_FILE"

echo "Thread study completed at $(date)" >> "$LOG_FILE"

##!/bin/bash
##SBATCH --nodes=1
##SBATCH --cpus-per-task=112
##SBATCH --ntasks-per-node=1
##SBATCH --mem=0
##SBATCH --time=00:30:00
##SBATCH --job-name=thread_scaling

#set -x
#DIR=$(pwd)/thread_exp_
#if [ ! -d "$DIR" ]; then
#  mkdir -p "$DIR"
#fi
## Initialize environment
#LOG_FILE="$DIR/thread_scaling.log"
#RESULTS_FILE="$DIR/thread_results.csv"
#SPEEDUP_FILE="$DIR/thread_speedup.csv"

## Set environment
#if [[ $ENV == "LEONARDO" ]]; then
#    echo "[$(date)] Setting environment for LEONARDO" >> "$LOG_FILE"
#    module load openmpi/4.1.6--gcc--12.2.0
#elif [[ $ENV == "ORFEO" ]]; then
#    echo "[$(date)] Setting environment for ORFEO" >> "$LOG_FILE"
#    module load openMPI
#fi

#echo "Started thread scalability study at $(date)" > "$LOG_FILE"

#EXEC=./stencil_parallel

## Initialize results file
#echo "threads,time" > "$RESULTS_FILE"

## Test different thread counts
#for threads in 1 2 4 8 16 32 56 84 112; do
#    echo "[$(date)] Testing OMP_NUM_THREADS=$threads" >> "$LOG_FILE"
    
#    export OMP_NUM_THREADS=$threads
#    export OMP_PLACES=cores
#    export OMP_PROC_BIND=close

#    # Run with timing
#    start_time=$(date +%s.%N)
#    /usr/bin/time -f "%e" -o time_tmp.txt $EXEC -o 1 > out_${threads}threads.txt 2> err_${threads}threads.txt
#    runtime=$(cat time_tmp.txt)
#    end_time=$(date +%s.%N)
#    wall_time=$(echo "$end_time - $start_time" | bc)
    
#    # Log results
#    echo "$threads,$runtime,$wall_time" >> "$RESULTS_FILE"
#    echo "Completed $threads threads in $runtime s (wall: $wall_time s)" >> "$LOG_FILE"
#done

## Compute speedup and efficiency
#awk -F, '
#  BEGIN {
#    print "threads,time,speedup,efficiency,wall_time"
#  }
#  NR==1 {next}  # Skip header
#  NR==2 {
#    t1 = $2  # Baseline time (1 thread)
#  }
#  {
#    speedup = t1/$2
#    efficiency = speedup/$1
#    printf("%d,%.2f,%.2f,%.2f,%.2f\n", $1, $2, speedup, efficiency, $3)
#  }
#' "$RESULTS_FILE" > "$SPEEDUP_FILE"

#echo "Thread study completed at $(date)" >> "$LOG_FILE"
