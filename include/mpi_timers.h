#ifndef MPI_TIMERS_H
#define MPI_TIMERS_H

#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <errno.h>

extern double comm_time;
extern double comp_time;
extern double loop_start_time;
extern double loop_end_time;
extern double factorization_time;

#define MAX_THREADS 128

#define NUM_TIMED_FUNCS 2 // compute and communication

extern double thread_times[NUM_TIMED_FUNCS][MAX_THREADS]; // Thread-local timings


#define MPI_CALL_TIMER(call, time_var)      \
    do {                                   \
        double t0 = MPI_Wtime();           \
        call;                              \
        time_var += MPI_Wtime() - t0;      \
    } while(0)


#define START_LOOP_TIMER() (loop_start_time = MPI_Wtime())
#define STOP_LOOP_TIMER()  (loop_end_time = MPI_Wtime())

// ---- Thread Timing Logger ---- //
static void log_thread_stats(int num_threads, int Rank, int func_id) {
    if (num_threads <= 1 || Rank != 0) return;

    double max_time = 0.0, min_time = 1e9, sum_time = 0.0;
    for (int t = 0; t < num_threads; ++t) {
        double tval = thread_times[func_id][t];
        if (tval > max_time) max_time = tval;
        if (tval < min_time) min_time = tval;
        sum_time += tval;
    }
    
    double imbalance_ratio = 100.0 * (max_time - min_time) / max_time;
    double avg_time = sum_time / num_threads;

    // Log
    printf("Func %d | Thread stats: Avg=%.6fs Max=%.6fs Min=%.6fs Imbalance=%.2f%%\n",
           func_id, avg_time, max_time, min_time, imbalance_ratio);
   char datetime_str[32];
   time_t now = time(NULL);
   struct tm *t = localtime(&now);
   strftime(datetime_str, sizeof(datetime_str), "%Y-%m-%d_%H%M", t);
   //
   // Parse node count from environment variable
   int num_nodes = 0;
   char *env_nodes = getenv("SLURM_NNODES");
   if (env_nodes != NULL){
	num_nodes = atoi(env_nodes);
    }

    char filename[64];    // allocate 

    // Log to CSV (appends if file exists)
    	
    snprintf(filename, sizeof(filename), "thread_stats_%d_%d.csv", num_threads, datetime_str);

    FILE *csv = fopen(filename, "a");
    if (csv) {
    // Write header if file is empty
    fseek(csv, 0, SEEK_END);
    if (ftell(csv) == 0) {  // File is empty
	fprintf(csv, "num_threads, num_nodes,func_id,avg_time,max_time,min_time,imbalance_ratio\n");
     }
     // Write data
    fprintf(csv, "%d, %d, %d,%.6f,%.6f,%.6f,%.2f\n",
	    num_threads, num_nodes, func_id, avg_time, max_time, min_time, imbalance_ratio);
    fclose(csv);
    } else {
	    fprintf(stderr, "Error opening file '%s': %s\n", filename, strerror(errno));
    }
    /* FILE *csv = fopen(filename, "a"); */
    /* fprintf("num_threads, func_id, avg_time, max_time, min_time, imbalance_ratio"); */
    /* if (csv) { */
    /*     fprintf(csv, "%d,%d,%.6f,%.6f,%.6f,%.2f\n", */
    /*             num_threads, func_id, avg_time, max_time, min_time, imbalance_ratio); */
    /*     fclose(csv); */
    /* } */
    /* else if (csv == NULL) { */
    /*     fprintf(stderr, "Error opening file '%s': %s\n", filename, strerror(errno)); */
    /*     perror("Error opening file"); */
    /* } */
}

// ---- Per-Rank Log Helper ---- //
static void log_rank_stats(int Rank, const char* label, double total_time, 
                          double comm_time, double compute_time, double mem_MB,
                          int num_threads) {
    char datetime_str[32];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    strftime(datetime_str, sizeof(datetime_str), "%Y-%m-%d_%H%M", t);

    char fname[64];
    snprintf(fname, sizeof(fname), "timing_rank_%d_%s.log", Rank, datetime_str);

    FILE *f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "Rank %d: Failed to open log file\n", Rank);
        return;
    }

    fprintf(f, "RANK %d STATS\n==============\n", Rank);
    fprintf(f, "Label:      %s\n", label);
    fprintf(f, "Total time: %.6f s\n", total_time);
    fprintf(f, "Comm time:  %.6f s\n", comm_time);
    fprintf(f, "Comp time:  %.6f s\n", compute_time);
    fprintf(f, "Mem usage:  %.2f MB\n", mem_MB);
    fprintf(f, "Threads:    %d\n\n", num_threads);

    for (int func_id = 0; func_id < NUM_TIMED_FUNCS; ++func_id) {
        fprintf(f, "Function %d thread timings:\n", func_id);
        for (int t = 0; t < num_threads; ++t) {
            fprintf(f, "  Thread %2d: %.6f s\n", t, thread_times[func_id][t]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

// ---- Main Reporting Function ---- //
static inline void report_timing_stats(MPI_Comm comm, int Rank, int Ntasks, 
                                     const char* label, int log_per_rank) {

    double total_time = loop_end_time - loop_start_time;
    double compute_time = total_time - comm_time;

    // Memory usage
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    double mem_MB = rusage.ru_maxrss / 1024.0;

    // Reductions
    double global_comm_time, max_comm_time, avg_comm_time;
    MPI_Reduce(&comm_time, &global_comm_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);//sum of times for comm
    MPI_Reduce(&comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);//max time for comm
    avg_comm_time = global_comm_time / Ntasks;//avg time for comm

    double global_comp_time, max_comp_time, min_comp_time;
    MPI_Reduce(&compute_time, &global_comp_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);//sum of times for comp
    MPI_Reduce(&compute_time, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);//max time for comp
    MPI_Reduce(&compute_time, &min_comp_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);//min time for comp
    double avg_comp_time = global_comp_time / Ntasks;//avg time for comp


    double imbalance_ratio = (max_comp_time - min_comp_time) / max_comp_time;

    double max_total_time;
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);//total time

    double total_factorization_time;
    MPI_Reduce(&factorization_time, &total_factorization_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);  

    int num_threads = omp_get_max_threads();
    
    if (Rank == 0) {
	// Generate datetime string for filename
	char datetime_str[32];
	time_t now = time(NULL);
	struct tm *t = localtime(&now);
	strftime(datetime_str, sizeof(datetime_str), "%Y-%m-%d_%H%M", t);

	char filename[128];
	snprintf(filename, sizeof(filename), "timing_report_%s.csv", datetime_str);
        


	FILE *csv = fopen(filename, "a");
	// Parse node count from environment variable
	int num_nodes = 0;
	char *env_nodes = getenv("SLURM_NNODES");

	if (env_nodes != NULL){
		num_nodes = atoi(env_nodes);
	}

        if (ftell(csv) == 0) {  // File is empty
		fprintf(csv, "label, num_nodes, threads,total_time,comm_time,avg_comm,avg_comp,max_comp,min_comp,imbalance_percent,tot_factorization_time\n");
        }

	
	if (csv) {
	fprintf(csv,
	    "%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f,%.6f\n",
	    label,
	    num_nodes,
	    num_threads,
	    max_total_time,
	    global_comm_time,
	    avg_comm_time,
	    avg_comp_time,
	    max_comp_time,
	    min_comp_time,
	    imbalance_ratio * 100,
	    total_factorization_time
	);
	fclose(csv);
	} else {
	perror("Could not open timing_report.csv for writing");
	}
        printf("\n=== %s Timing Report ===\n", label);
        printf("Wall clock time (max over ranks): %.6f s\n", max_total_time);
        printf("Total COMM time: %.6f s\n", global_comm_time);
        printf("COMP imbalance: %.2f%%\n", imbalance_ratio * 100);
	printf("total time spent for factorization with option: %.6f s\n", total_factorization_time);
    }

    // Thread stats (only logged by rank 0)
    for (int func_id = 0; func_id < NUM_TIMED_FUNCS; ++func_id) {
        log_thread_stats(num_threads, Rank, func_id);
    }

    // Per-rank logs (optional)
    if (log_per_rank) {
        log_rank_stats(Rank, label, total_time, comm_time, 
                      compute_time, mem_MB, num_threads);
    }
}
#endif
