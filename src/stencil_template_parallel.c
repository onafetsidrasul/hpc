/*

/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */


#include "stencil_template_parallel.h"

double comm_time = 0.0, comp_time = 0.0;
double loop_start_time = 0.0, loop_end_time = 0.0;
double thread_times[NUM_TIMED_FUNCS][MAX_THREADS] = {1.0};
double factorization_time = 0.0;

// ------------------------------------------------------------------
// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  MPI_Comm myCOMM_WORLD;
  int  Rank, Ntasks;
  uint neighbours[4];
  int verbose = 0;

  int  Niterations;
  int  periodic;
  vec2_t S, N;
  
  int      Nsources;
  int      Nsources_local;
  vec2_t  *Sources_local;
  double   energy_per_source;

  plane_t   planes[2];  
  buffers_t buffers[2];
  
  int output_energy_stat_perstep;
  
  /* initialize MPI envrionment */
  {
    int level_obtained;
    
    // NOTE: change MPI_FUNNELED if appropriate
    //
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
    if ( level_obtained < MPI_THREAD_FUNNELED ) {
      printf("MPI_thread level obtained is %d instead of %d\n",
	     level_obtained, MPI_THREAD_FUNNELED );
      MPI_Finalize();
      exit(1); }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);
  }
  
  
  /* argument checking and setting */
  int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
			 neighbours, &Niterations,
			 &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
			 &planes[0], &buffers[0], &verbose);

  if ( ret )
    {
      printf("task %d is opting out with termination code %d\n",
	     Rank, ret );
       fflush(stdout);
      
      MPI_Finalize();
      return 0;
    }
  
  
  int current = OLD;
  START_LOOP_TIMER();
  for (int iter = 0; iter < Niterations; ++iter)
    
    {
      
      
      
      /* new energy from sources */
      int inj_energy=inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, N, &planes[current], verbose);
      if ( inj_energy != 0 ) {
        fprintf(stderr, "Rank %d: error %d in inject_energy\n", Rank, inj_energy);
        MPI_Abort(MPI_COMM_WORLD, inj_energy);
      }
      else {
        if (verbose > 0)
          printf("Rank %d: inject_energy success\n", Rank);
          fflush(stdout);
      }

      /* -------------------------------------- */

      // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position
      const int register width = planes[current].size[_x_];
      const int register height = planes[current].size[_y_];
      const int register buffer_width = width + 2;
      const int register buffer_height = height + 2;

      pack_boundary( &planes[current], buffers, width, height, neighbours, verbose, Rank );
      // [B] perfoem the halo communications
      //     (1) use Send / Recv
      //     (2) use Isend / Irecv
      //         --> can you overlap communication and compution in this way?
      send_boundary(buffers, neighbours, buffer_width, buffer_height, Rank, verbose, 1); 
      
      // [C] copy the haloes data
      update_boundary( &planes[current], buffers, width, height, neighbours, verbose, Rank );
      /* --------------------------------------  */
      /* update grid points */
      
      update_plane( periodic, N, &planes[current], &planes[!current] );

      /* output if needed */
      if ( output_energy_stat_perstep )
	output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
      /* swap plane indexes for the new iteration */
      current = !current;
      
    }

  STOP_LOOP_TIMER();

  output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  if (Ntasks > 1){
  	MPI_CALL_TIMER(MPI_Barrier(myCOMM_WORLD), comm_time);
  }

  memory_release(buffers, planes, Rank, verbose);

  report_timing_stats(myCOMM_WORLD, Rank, Ntasks, (Ntasks == 1) ? "Serial Mode" : "Main Loop", 0);
  
  MPI_Finalize();
  return 0;
}