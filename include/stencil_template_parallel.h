/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 * See COPYRIGHT in top-level directory.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <math.h>

#include <omp.h>
#include <mpi.h>

#include "mpi_timers.h"

#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1

#define TAG_BORDER_EXCHANGE 77 

// scheduling types for OpenMP
#define STATIC   0 // Good for loops with similar workloads per iteration.
#define DYNAMIC  1 // Good for loops with varying workloads per iteration.
#define GUIDED   2 // Good for loops with unpredictable workloads.
#define AUTO     3 // Let OpenMP decide the best scheduling.

#ifndef OPENMP_SCHEDULE
#define OPENMP_SCHEDULE STATIC
#endif

typedef unsigned int uint;

typedef uint    vec2_t[2];
typedef double *restrict buffers_t[4];

typedef struct {
    double   * restrict data;
    vec2_t     size;
} plane_t;




extern int inject_energy ( const int      ,
                           const int      ,
			                const vec2_t  *,
			                const double   ,
                           const vec2_t   ,   
                           plane_t *      ,  
			            int             );


extern int update_plane ( const int      ,
                          const vec2_t   ,
                          const plane_t *,
                                plane_t * );


extern int get_total_energy( plane_t *,
                             double  * );


// for the domain decomposition!! Trying to find an optimal one with less communication surface
// Less messages, less latency
// It is serial but it's fine
uint simple_factorization(uint A, int *Nfactors, uint **factors) {
    int max_threads = omp_get_max_threads();
    uint **local_factors = malloc(max_threads * sizeof(uint *));
    int *local_counts = calloc(max_threads, sizeof(int));
    int *local_capacities = malloc(max_threads * sizeof(int));

    for (int i = 0; i < max_threads; ++i) {
        local_capacities[i] = 16;
        local_factors[i] = malloc(local_capacities[i] * sizeof(uint));
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int local_count = 0;
        uint local_A = A;

        #if OPENMP_SCHEDULE == SCHED_STATIC
            #pragma omp for schedule(static, OMP_CHUNK_SIZE) nowait
        #elif OPENMP_SCHEDULE == SCHED_DYNAMIC
            #pragma omp for schedule(dynamic, OMP_CHUNK_SIZE) nowait
        #elif OPENMP_SCHEDULE == SCHED_GUIDED
            #pragma omp for schedule(guided, OMP_CHUNK_SIZE) nowait
        #else
            #pragma omp for schedule(auto) nowait
        #endif
        for (int f = 2; f <= A / 2; ++f) {
            while (local_A % f == 0) {
                if (local_count >= local_capacities[tid]) {
                    local_capacities[tid] *= 2;
                    local_factors[tid] = realloc(local_factors[tid], local_capacities[tid] * sizeof(uint));
                    if (!local_factors[tid]) {
                        fprintf(stderr, "Error: realloc failed\n");
                        exit(EXIT_FAILURE);
                    }
                }
                local_factors[tid][local_count++] = f;
                local_A /= f;
            }
        }

        local_counts[tid] = local_count;
    }

    int total_count = 0;
    for (int i = 0; i < max_threads; ++i)
        total_count += local_counts[i];

    uint *final_factors = malloc(total_count * sizeof(uint));
    int offset = 0;
    for (int i = 0; i < max_threads; ++i) {
        for (int j = 0; j < local_counts[i]; ++j)
            final_factors[offset++] = local_factors[i][j];
        free(local_factors[i]);
    }

    free(local_factors);
    free(local_counts);
    free(local_capacities);

    *Nfactors = total_count;
    *factors = final_factors;
    return 0;
}


extern int memory_allocate ( const int       *,
 		      const vec2_t     ,
 		      buffers_t *      ,
 		      plane_t   *      );


extern int memory_release (buffers_t *, plane_t   *,int ,  int );


int output_energy_stat ( int      ,
                         plane_t *,
                         double   ,
                         int      ,
                         MPI_Comm *);




inline int inject_energy ( const int      periodic,
                           const int      Nsources,
			               const vec2_t  *Sources,
			               const double   energy,
                           const vec2_t   N,
                           plane_t *plane,
                           int verbose)
{
    const uint register sizex = plane->size[_x_] + 2;
    const uint register sizey = plane->size[_y_];
    if (sizex <= 0 || sizey <= 0) {
        fprintf(stderr, "ERR: Invalid size: x=%d, y=%d\n", plane->size[_x_], plane->size[_y_]);
        exit(EXIT_FAILURE);
    }
    if (verbose > 0){
	    printf("sizex: %d sizey: %d \n", sizex, sizey);
	    printf("total data size (sizex * sizey %d \n", sizex*sizey);
	    fflush(stdout);
    }

    // only pointer that access the data, it allows better optimization
    double * restrict data = plane->data;

    #define IDX( i, j ) ( (j)*sizex + 2 + (i) )
    for (int s = 0; s < Nsources; s++)
        {
            int x = Sources[s][_x_];
            int y = Sources[s][_y_];
	    if (verbose > 0){
	    
		    printf("Source %d: x = %d, y = %d\n, Nsources: %d \n", s, x, y, Nsources);
		    fflush(stdout);

		    printf("IDX(%d, %d) = %d \n", x, y, IDX(x, y));
		    fflush(stdout);
	    }

            data[ IDX(x,y) ] += energy;
            
            if ( periodic )
                {
                    if ( (N[_x_] == 1)  )
                        {
                            // propagate the boundaries 
                            data[IDX(N[_x_]+1, y)] += energy;
                        }
                    
                    if ( (N[_y_] == 1) )
                        {
                            // propagate the boundaries
                            data[IDX(x, N[_y_]+1)] += energy;
                        }
                }                
        }
 #undef IDX
    
  return 0;
}


inline int update_plane ( const int      periodic, 
                          const vec2_t   N,         // the grid of MPI tasks
                          const plane_t *oldplane,
                                plane_t *newplane
                          )
    
{
    uint fxsize = oldplane->size[_x_] + 2;
    uint fysize = oldplane->size[_y_] + 2;
    uint xsize = oldplane->size[_x_];
    uint ysize = oldplane->size[_y_];

    #define IDX(i, j) ((j) * fxsize + (i) )

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    double reduced_sum = 0.0;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        double t_start = omp_get_wtime();

        double summ = 0.0;

        // Select schedule policy via preprocessor (I can choose it)
        // nowait to avoid implicit barrier, the thread can proceed without waiting for the others
        #if OPENMP_SCHEDULE == STATIC
            #pragma omp for schedule(static) nowait
        #elif OPENMP_SCHEDULE == DYNAMIC
            #pragma omp for schedule(dynamic, OMP_CHUNK_SIZE) nowait
        #elif OPENMP_SCHEDULE == GUIDED
            #pragma omp for schedule(guided, OMP_CHUNK_SIZE) nowait
        #else //Auto
            #pragma omp for schedule(auto) nowait
        #endif

        for (uint j = 1; j <= ysize; j++) {
            for (uint i = 1; i <= xsize; i++) {
                double alpha = 0.6;
                double result = old[IDX(i,j)] * alpha;
                double sum_i = (old[IDX(i-1, j)] + old[IDX(i+1, j)]) / 4.0 * (1 - alpha);
                double sum_j = (old[IDX(i, j-1)] + old[IDX(i, j+1)]) / 4.0 * (1 - alpha);
                result += (sum_i + sum_j);
                new[IDX(i,j)] = result;
                summ += result;
            }
        }

        #pragma omp atomic
        reduced_sum += summ;

        double t_end = omp_get_wtime();
        thread_times[1][tid] += (t_end - t_start);
    }

    if (periodic) {
        if (N[_x_] == 1) {
            for (int i = 1; i <= xsize; i++) {
                new[i] = new[IDX(i, ysize)];
                new[IDX(i, ysize+1)] = new[i];
            }
        }
        if (N[_y_] == 1) {
            for (int j = 1; j <= ysize; j++) {
                new[IDX(0, j)] = new[IDX(xsize, j)];
                new[IDX(xsize+1, j)] = new[IDX(1, j)];
            }
        }
    }

    return 0;
    #undef IDX
}




inline int get_total_energy( plane_t *plane,
                             double  *energy )
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{

    const int register xsize = plane->size[_x_];
    const int register ysize = plane->size[_y_];

    const int register fsize = xsize+2;

    double * restrict data = plane->data;
    
   #define IDX( i, j ) ( (j)*fsize + (i) )

   double totenergy = 0;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        double start = omp_get_wtime();

        double private_energy = 0.0;
        #if OPENMP_SCHEDULE == STATIC
            #pragma omp for schedule(static, OMP_CHUNK_SIZE) nowait
        #elif OPENMP_SCHEDULE == DYNAMIC
            #pragma omp for schedule(dynamic, OMP_CHUNK_SIZE) nowait
        #elif OPENMP_SCHEDULE == GUIDED
            #pragma omp for schedule(guided, OMP_CHUNK_SIZE) nowait
        #else  // auto
            #pragma omp for schedule(auto) nowait
        #endif
        
        for (int j = 1; j <= ysize; j++) {
            for (int i = 1; i <= xsize; i++) {
                private_energy += data[IDX(i, j)];
            }
        }

        #pragma omp atomic
        totenergy += private_energy;

        double end = omp_get_wtime();
        thread_times[0][tid] += (end - start);
    }
    *energy = totenergy;
    return 0;
    #undef IDX
}

inline int memory_allocate ( const int       *neighbours  ,
		            const vec2_t     N           ,
		            buffers_t *buffers_ptr ,
		            plane_t   *planes_ptr
		            )

{
    /*
      here you allocate the memory buffers that you need to
      (i)  hold the results of your computation
      (ii) communicate with your neighbours

      The memory layout that I propose to you is as follows:

      (i) --- calculations
      you need 2 memory regions: the "OLD" one that contains the
      results for the step (i-1)th, and the "NEW" one that will contain
      the updated results from the step ith.

      Then, the "NEW" will be treated as "OLD" and viceversa.

      These two memory regions are indexed by *plane_ptr:

      planes_ptr[0] ==> the "OLD" region
      plames_ptr[1] ==> the "NEW" region
      


      (ii) --- communications

      you may need two buffers (one for sending and one for receiving)
      for each one of your neighnours, that are at most 4:
      north, south, east amd west.      

      To them you need to communicate at most mysizex or mysizey
      double data.

      These buffers are indexed by the buffer_ptr pointer so
      that

      (*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
      (*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
      
      --->> Of course you can change this layout as you prefer
      
     */

  if (planes_ptr == NULL )
      {
       int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            fprintf(stderr, "[Rank %d] Fatal: NULL planes pointer\n", rank);
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      } 
 

  if (buffers_ptr == NULL )
    {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        fprintf(stderr, "[Rank %d] Fatal: NULL buffers pointer\n", rank);
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    

  // ··················································
  // allocate memory for data
  // we allocate the space needed for the plane plus a contour frame
  // that will contain data form neighbouring MPI tasks

  unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2);

  planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[OLD].data == NULL ){
        fprintf(stderr, "Error: Memory allocation for planes_ptr failed\n");
        exit(EXIT_FAILURE); // or return an error code
  }
  memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) );

  planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[NEW].data == NULL ){
    // manage the malloc fail
        fprintf(stderr, "Error: Memory allocation for planes_ptr failed\n");
        exit(EXIT_FAILURE); // or return an error code
  }
  memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );


  // ··················································
  // buffers for north and south communication 
  // are not really needed
  //
  // in fact, they are already contiguous, just the
  // first and last line of every rank's plane
  //
  // you may just make some pointers pointing to the
  // correct positions
  //

  // or, if you preer, just go on and allocate buffers
  // also for north and south communications

  // ··················································
  // allocate buffers
  // unsigned int buffer_frame_size = (buffers_ptr[OLD].size[_x_]+2) * (buffers_ptr[OLD].size[_y_]+2);

  unsigned int buffer_frame_size = frame_size; 

  // allocate buffers 
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 2; ++j){
    (buffers_ptr)[j][i] = malloc(frame_size * sizeof(double));
    if ((buffers_ptr)[j][i] == NULL) {
        fprintf(stderr, "Error: Memory allocation for (*buffers_ptr)[%d][%d] failed\n", j, i);
        exit(EXIT_FAILURE);
    }
    memset((buffers_ptr)[j][i], 0, frame_size * sizeof(double));
    }

  // ··················································
  
  return 0;
}

inline int memory_release(buffers_t *buffers, plane_t *planes, int Rank, int verbose) {

    if (planes[OLD].data) {
        free(planes[OLD].data);
        planes[OLD].data = NULL;
    }
    if (planes[NEW].data) {
        free(planes[NEW].data);
        planes[NEW].data = NULL;
    }
   // free the communication buffers
    if (buffers) {
        for (int i = 0; i < 4; i++) {
            buffers[i][0] = NULL;
        
	    if (Rank == 0)
		    printf("plane at rank %d freed.\n", Rank);

	}
    }
    return 0;
}

extern int initialize_sources( int       ,
			int       ,
			MPI_Comm  *,
			uint      [2],
			int       ,
			int      *,
			vec2_t  ** );

inline int initialize_sources( int       Me,
			int       Ntasks,
			MPI_Comm *Comm,
			vec2_t    mysize,
			int       Nsources,
			int      *Nsources_local,
			vec2_t  **Sources )

{

  srand48(time(NULL) ^ Me);
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  // Master randomly generates an assignment of sources to tasks
  if ( Me == 0 )
    {
      for ( int i = 0; i < Nsources; i++ )
	tasks_with_sources[i] = (int)lrand48() % Ntasks;
    }
 
  // Broadcast the assignments to all tasks
  MPI_CALL_TIMER(MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm ), comm_time);

  // How many sources are local to this task?
  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ )
    nlocal += (tasks_with_sources[i] == Me);
  *Nsources_local = nlocal;
  
  // Allocate and fill the local sources array
  if ( nlocal > 0 )
    {
      vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
      for ( int s = 0; s < nlocal; s++ )
	{
	  helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	  helper[s][_y_] = 1 + lrand48() % mysize[_y_];
	}

      *Sources = helper;
    }
  
  free( tasks_with_sources );

  return 0;
}

extern int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm );

inline int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm )
{

  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  // aggregate the total energy
  MPI_CALL_TIMER(MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm ), comm_time);

  // Print the energy statistics (only on rank 0)
  if ( Me == 0 )
    {
      if ( step >= 0 )
	printf(" [ step %4d ] ", step ); fflush(stdout);

      
      printf( "total injected energy is %g, "
	      "system energy is %g "
	      "( in avg %g per grid point)\n",
	      budget,
	      tot_system_energy,
	      tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
    }
  
  return 0;
}


extern int initialize ( MPI_Comm *,
                 int       ,
		 int       ,
		 int       ,
		 char    **,
                 vec2_t   *,
                 vec2_t   *,                 
		 int      *,
                 int      *,
		 int      *,
		 int      *,
		 int      *,
		 int      *,
                 vec2_t  **,
                 double   *,
                 plane_t  *,
                 buffers_t *,
		 int       *);

inline int initialize ( MPI_Comm *Comm,
		 int      Me,                  // the rank of the calling process
		 int      Ntasks,              // the total number of MPI ranks
		 int      argc,                // the argc from command line
		 char   **argv,                // the argv from command line
		 vec2_t  *S,                   // the size of the plane
		 vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
		 int     *periodic,            // periodic-boundary tag
		 int     *output_energy_stat,
		 int     *neighbours,          // four-int array that gives back the neighbours of the calling task
		 int     *Niterations,         // how many iterations
		 int     *Nsources,            // how many heat sources
		 int     *Nsources_local,
		 vec2_t **Sources_local,
		 double  *energy_per_source,   // how much heat per source
		 plane_t *planes,
		 buffers_t *buffers,
		 int      * verbose)
{
  int halt = 0;
  int ret;
  
  // ··································································
  // set deffault values

  (*S)[_x_]         = 10000;
  (*S)[_y_]         = 10000;
  *periodic         = 0;
  *Nsources         = 4;
  *Nsources_local   = 0;
  *Sources_local    = NULL;
  *Niterations      = 1000;
  *energy_per_source = 1.0;
  *verbose = 0;
  if ( planes == NULL ) {
    // manage the situation
  }

  planes[OLD].size[0] = planes[OLD].size[0] = 0;
  planes[NEW].size[0] = planes[NEW].size[0] = 0;
  
  for ( int i = 0; i < 4; i++ )
    neighbours[i] = MPI_PROC_NULL;

  for ( int b = 0; b < 2; b++ )
    for ( int d = 0; d < 4; d++ )
      buffers[b][d] = NULL;
  
  // ··································································
  // process the commadn line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":h:x:y:e:E:n:o:p:v:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': (*S)[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': (*S)[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'o': *output_energy_stat = (atoi(optarg) > 0);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'v': *verbose = atoi(optarg);
	    break;

	  case 'h': {
	    if ( Me == 0 )
	      printf( "\nvalid options are ( values btw [] are the default values ):\n"
		      "-x    x size of the plate [10000]\n"
		      "-y    y size of the plate [10000]\n"
		      "-e    how many energy sources on the plate [4]\n"
		      "-E    how many energy sources on the plate [1.0]\n"
		      "-n    how many iterations [1000]\n"
		      "-p    whether periodic boundaries applies  [0 = false]\n\n"
		      );
	    halt = 1; }
	    break;
	    
	    
	  case ':': printf( "option -%c requires an argument\n", optopt);
	    break;
	    
	  case '?': printf(" -------- help unavailable ----------\n");
	    break;
	  }
      }

    if ( opt == -1 )
      break;
  }

  if ( halt )
    return 1;
  
  
  // ··································································
  /*
   * here we should check for all the parms being meaningful
   *
   */
// Validate parameters
int valid = 1;
if ((*S)[_x_] == 0 || (*S)[_y_] == 0) {
    fprintf(stderr, "Rank %d: Error - Plate dimensions cannot be zero\n", Me);
    valid = 0;
}

if (*Nsources < 0) {
    fprintf(stderr, "Rank %d: Error - Number of sources cannot be negative\n", Me);
    valid = 0;
}

if (*Niterations <= 0) {
    fprintf(stderr, "Rank %d: Error - Number of iterations must be positive\n", Me);
    valid = 0;
}

if (*energy_per_source <= 0) {
    fprintf(stderr, "Rank %d: Error - Energy per source must be positive\n", Me);
    valid = 0;
}


if (!valid) {
    MPI_Abort(*Comm, 1);
    return 1;
}

if (Me == 0 && *verbose > 0) {
    printf("\nParameter validation and values:\n");
    printf("--------------------------------\n");
    printf("Plate size (x, y): %u, %u\n", (*S)[_x_], (*S)[_y_]);
    printf("Periodic boundaries: %d\n", *periodic);
    printf("Number of heat sources: %d\n", *Nsources);
    printf("Energy per source: %.2f\n", *energy_per_source);
    printf("Number of iterations: %d\n", *Niterations);
    printf("Output energy statistics: %d\n", *output_energy_stat);
    printf("Verbose level: %d\n", *verbose);
    printf("--------------------------------\n\n");
    fflush(stdout);
}
  // ...

  
  // ··································································
  /*
   * find a suitable domain decomposition
   * very simple algorithm, you may want to
   * substitute it with a better one
   *
   * the plane Sx x Sy will be solved with a grid
   * of Nx x Ny MPI tasks
   */

  vec2_t Grid;
  double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
  int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) );

  // if the number of tasks is less than the formfactor + 1
  // we can use only one dimension for the decomposition
  
  if ( dimensions == 1 )
    {
      if ( (*S)[_x_] >= (*S)[_y_] )
	Grid[_x_] = Ntasks, Grid[_y_] = 1;
      else
	Grid[_x_] = 1, Grid[_y_] = Ntasks;
    }
  else
    {
      int   Nf;
      uint *factors;
      uint  first = 1;
      double start_fact = omp_get_wtime();     
      ret = simple_factorization( Ntasks, &Nf, &factors );
      double end_fact = omp_get_wtime();
      factorization_time = end_fact - start_fact;  
      for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ )
	first *= factors[i];

      if ( (*S)[_x_] > (*S)[_y_] )
	Grid[_x_] = Ntasks/first, Grid[_y_] = first;
      else
	Grid[_x_] = first, Grid[_y_] = Ntasks/first;
    }

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];
  

  // ··································································
  // my cooridnates in the grid of processors
  //
  int X = Me % Grid[_x_];
  int Y = Me / Grid[_x_];

  // ··································································
  // find my neighbours
  //

  if ( Grid[_x_] > 1 )
    {  
      if ( *periodic ) {       
	neighbours[EAST]  = Y*Grid[_x_] + (Me + 1 ) % Grid[_x_];
	neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
      
      else {
	neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
	neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
    }

  if ( Grid[_y_] > 1 )
    {
      if ( *periodic ) {      
	neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; }

      else {    
	neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
	neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
    }

  // ··································································
  // the size of my patch
  //

  /*
   * every MPI task determines the size sx x sy of its own domain
   * REMIND: the computational domain will be embedded into a frame
   *         that is (sx+2) x (sy+2)
   *         the outern frame will be used for halo communication or
   */
  
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_];
  uint r = (*S)[_x_] % Grid[_x_];
  mysize[_x_] = s + (X < r);
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r);

  planes[OLD].size[0] = mysize[0];
  planes[OLD].size[1] = mysize[1];
  planes[NEW].size[0] = mysize[0];
  planes[NEW].size[1] = mysize[1];
  

  if ( *verbose > 0 )
    {
      if ( Me == 0 ) {
	printf("Tasks are decomposed in a grid %d x %d\n\n",
		 Grid[_x_], Grid[_y_] );
	fflush(stdout);
      }

      MPI_CALL_TIMER(MPI_Barrier(*Comm), comm_time);
      
      for ( int t = 0; t < Ntasks; t++ )
	{
	  if ( t == Me )
	    {
	      printf("Task %4d :: "
		     "\tgrid coordinates : %3d, %3d\n"
		     "\tneighbours: N %4d    E %4d    S %4d    W %4d\n",
		     Me, X, Y,
		     neighbours[NORTH], neighbours[EAST],
		     neighbours[SOUTH], neighbours[WEST] );
	      fflush(stdout);
	    }

	  MPI_Barrier(*Comm);
	}
      printf("neighbours determined and barrier overcome for Task %d.\n Onto Memory initialization now \n", Me);
      fflush(stdout);
    }

  
  // ··································································
  // allocate the needed memory
  //
  ret = memory_allocate(neighbours,
                       *N, 
                       buffers, 
                       planes); 
  if (ret != 0) {
       fprintf(stderr, "Rank %d: memory_allocate failed with code %d\n", Me, ret);
       MPI_Abort(MPI_COMM_WORLD, ret);  
  } else {
        printf("Rank %d: memory_allocate succeeded\n", Me);
        fflush(stdout);
	}	
  // ··································································
  // allocate the heat sources
  //
  ret = initialize_sources( Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local );

  if (ret != 0) {
    fprintf(stderr, "Rank %d: initialize_sources failed with code %d\n", Me, ret);
    MPI_Abort(MPI_COMM_WORLD, ret);
  } else {
    printf("Rank %d: initialize_sources succeeded\n", Me);
    fflush(stdout);
  }
  return 0;  
}

/**
 *
 * Updates ghost cells in the plane using received halo data (buffers[RECV][DIR]).
 * 
 * @param plane          The plane to update (planes[current] or planes[!current]).
 * @param buffers        Array of send/recv buffers (buffers[RECV][DIR] for unpacking).
 * @param width         Width of the local subdomain (including halos).
 * @param height        Height of the local subdomain (including halos).
 * @param neighbours    Array of neighbour ranks (MPI_PROC_NULL if no neighbour).
 */
void update_boundary(plane_t *plane, buffers_t buffers[2], 
                      int width, int height, const int neighbours[4],
                      int verbose, int rank) {
    
    if (verbose > 0) {
        printf("Rank %d: Updating ghost cells...\n", rank);
    }

    // ---- EAST (left ghost column) ----
    if (neighbours[EAST] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Updating EAST ghost cells (left column at x=0):\n", rank);
        }
        for (int j = 1; j < height - 1; j++) {
            if (verbose > 0) {
                printf("  y=%d: old=%f, new=%f\n", 
                       j, plane->data[j * width], buffers[RECV][EAST][j - 1]);
            }
            plane->data[j * width] = buffers[RECV][EAST][j - 1];
        }
    }

    // ---- WEST (right ghost column) ----
    if (neighbours[WEST] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Updating WEST ghost cells (right column at x=%d):\n", rank, width - 1);
        }
        for (int j = 1; j < height - 1; j++) {
            if (verbose > 0) {
                printf("  y=%d: old=%f, new=%f\n", 
                       j, plane->data[(j + 1) * width - 1], buffers[RECV][WEST][j - 1]);
            }
            plane->data[(j + 1) * width - 1] = buffers[RECV][WEST][j - 1];
        }
    }

    // ---- NORTH (top ghost row) ----
    if (neighbours[NORTH] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Updating NORTH ghost cells (top row at y=0):\n", rank);
        }
        for (int i = 1; i < width - 1; i++) {
            if (verbose > 0) {
                printf("  x=%d: old=%f, new=%f\n", 
                       i, plane->data[i], buffers[RECV][NORTH][i - 1]);
            }
            plane->data[i] = buffers[RECV][NORTH][i - 1];
        }
    }

    // ---- SOUTH (bottom ghost row) ----
    if (neighbours[SOUTH] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Updating SOUTH ghost cells (bottom row at y=%d):\n", rank, height - 1);
        }
        for (int i = 1; i < width - 1; i++) {
            if (verbose > 0) {
                printf("  x=%d: old=%f, new=%f\n", 
                       i, plane->data[(height - 1) * width + i], buffers[RECV][SOUTH][i - 1]);
            }
            plane->data[(height - 1) * width + i] = buffers[RECV][SOUTH][i - 1];
        }
    }
}

/**
 * Packs boundary data into send buffers (buffers[SEND][DIR]).
 * 
 * @param plane          The plane to read from (planes[current]).
 * @param buffers        Array of send/recv buffers (buffers[SEND][DIR] for packing).
 * @param width         Width of the local subdomain (including halos).
 * @param height        Height of the local subdomain (including halos).
 * @param neighbours    Array of neighbour ranks (MPI_PROC_NULL if no neighbour).
 */

void pack_boundary(const plane_t *plane, buffers_t buffers[2], 
                int width, int height, const int neighbours[4],
                int verbose, int rank) {
    
    if (verbose > 0) {
        printf("Rank %d: Packing halos...\n", rank);
    	fflush(stdout);

    }

    // ---- EAST (send left interior column) ----
    if (neighbours[EAST] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Packing EAST halo (interior column at x=1):\n", rank);
    	fflush(stdout);
        }
        for (int j = 1; j < height - 1; j++) {
            buffers[SEND][EAST][j - 1] = plane->data[j * width + 1];
            if (verbose > 0) {
                printf("  y=%d: %f\n", j, plane->data[j * width + 1]);
    		fflush(stdout);
            }
        }
    }

    // ---- WEST (send right interior column) ----
    if (neighbours[WEST] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Packing WEST halo (interior column at x=%d):\n", rank, width - 2);
    	fflush(stdout);
        }
        for (int j = 1; j < height - 1; j++) {
            buffers[SEND][WEST][j - 1] = plane->data[(j + 1) * width - 2];
            if (verbose> 0) {
                printf("  y=%d: %f\n", j, plane->data[(j + 1) * width - 2]);
    		fflush(stdout);
            }
        }
    }

    // ---- NORTH (send top interior row) ----
    if (neighbours[NORTH] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Packing NORTH halo (interior row at y=1):\n", rank);
    	fflush(stdout);
        }
        for (int i = 1; i < width - 1; i++) {
            buffers[SEND][NORTH][i - 1] = plane->data[width + i];
            if (verbose> 0) {
                printf("  x=%d: %f\n", i, plane->data[width + i]);
    		fflush(stdout);
            }
        }
    }

    // ---- SOUTH (send bottom interior row) ----
    if (neighbours[SOUTH] != MPI_PROC_NULL) {
        if (verbose > 0) {
            printf("Rank %d: Packing SOUTH halo (interior row at y=%d):\n", rank, height - 2);
    	fflush(stdout);
        }
        for (int i = 1; i < width - 1; i++) {
            buffers[SEND][SOUTH][i - 1] = plane->data[(height - 2) * width + i];
            if (verbose > 0) {
                printf("  x=%d: %f\n", i, plane->data[(height - 2) * width + i]);
    		fflush(stdout);
            }
        }
    }
}


void send_boundary(buffers_t buffers[2], const int neighbours[4], int buffer_width, int buffer_height, int Rank, int verbose, int non_blocking) {
    if (non_blocking) {
        MPI_Request reqs[8];
        int req_idx = 0;

        // EAST-WEST Comm
        if (neighbours[EAST] != MPI_PROC_NULL) {
            if (verbose > 0) {
                printf("Rank %d: Sending EAST  to %d:elements %d.(bytes: %d) \n", Rank, neighbours[EAST], buffer_height, buffer_height*sizeof(double));
            }
            MPI_CALL_TIMER(MPI_Isend(buffers[SEND][EAST], buffer_height, MPI_DOUBLE, neighbours[EAST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }
        if (neighbours[WEST] != MPI_PROC_NULL) {
	    if (verbose > 0) {
		printf("Rank %d: Receiving from WEST (%d) — expecting %d elements (bytes: %zu)\n",
		       Rank, neighbours[WEST], buffer_height, buffer_height * sizeof(double));
	    }
            MPI_CALL_TIMER(MPI_Irecv(buffers[RECV][WEST], buffer_height, MPI_DOUBLE, neighbours[WEST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }

        if (neighbours[WEST] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Isend(buffers[SEND][WEST], buffer_height, MPI_DOUBLE, neighbours[WEST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }
        if (neighbours[EAST] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Irecv(buffers[RECV][EAST], buffer_height, MPI_DOUBLE, neighbours[EAST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }

        // NORTH-SOUTH Comm
        if (neighbours[NORTH] != MPI_PROC_NULL) {
            if (verbose > 0) {
                printf("Rank %d: Sending NORTH to %d.\n", Rank, neighbours[NORTH]);
            }
            MPI_CALL_TIMER(MPI_Isend(buffers[SEND][NORTH], buffer_width, MPI_DOUBLE, neighbours[NORTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }
        if (neighbours[SOUTH] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Irecv(buffers[RECV][SOUTH], buffer_width, MPI_DOUBLE, neighbours[SOUTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }

        if (neighbours[SOUTH] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Isend(buffers[SEND][SOUTH], buffer_width, MPI_DOUBLE, neighbours[SOUTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }
        if (neighbours[NORTH] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Irecv(buffers[RECV][NORTH], buffer_width, MPI_DOUBLE, neighbours[NORTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, &reqs[req_idx++]), comm_time);
        }

        MPI_Waitall(req_idx, reqs, MPI_STATUSES_IGNORE);
    } else {
        // Blocking mode
        if (neighbours[EAST] != MPI_PROC_NULL) {
            if (verbose > 0) {
                printf("Rank %d: Sending EAST to %d.\n", Rank, neighbours[EAST]);
            }
            MPI_CALL_TIMER(MPI_Send(buffers[SEND][EAST], buffer_height, MPI_DOUBLE, neighbours[EAST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD), comm_time);
        }
        if (neighbours[WEST] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Recv(buffers[RECV][WEST], buffer_height, MPI_DOUBLE, neighbours[WEST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE), comm_time);
        }

        if (neighbours[WEST] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Send(buffers[SEND][WEST], buffer_height, MPI_DOUBLE, neighbours[WEST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD), comm_time);
        }
        if (neighbours[EAST] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Recv(buffers[RECV][EAST], buffer_height, MPI_DOUBLE, neighbours[EAST], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE), comm_time);
        }

        if (neighbours[NORTH] != MPI_PROC_NULL) {
            if (verbose > 0) {
                printf("Rank %d: Sending NORTH to %d.\n", Rank, neighbours[NORTH]);
            }
            MPI_CALL_TIMER(MPI_Send(buffers[SEND][NORTH], buffer_width, MPI_DOUBLE, neighbours[NORTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD), comm_time);
        }
        if (neighbours[SOUTH] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Recv(buffers[RECV][SOUTH], buffer_width, MPI_DOUBLE, neighbours[SOUTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE), comm_time);
        }

        if (neighbours[SOUTH] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Send(buffers[SEND][SOUTH], buffer_width, MPI_DOUBLE, neighbours[SOUTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD), comm_time);
        }
        if (neighbours[NORTH] != MPI_PROC_NULL) {
            MPI_CALL_TIMER(MPI_Recv(buffers[RECV][NORTH], buffer_width, MPI_DOUBLE, neighbours[NORTH], TAG_BORDER_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE), comm_time);
        }
    }
}