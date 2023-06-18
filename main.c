#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "variables.h"
#include "cosmoparam.h"
#include "correlacion.h"
#include "leefile.h"
#include "timer.h"

int main(int argc, char **argv)
{
  int NumProc, MyProc;
  double start_tot,end_tot,start,end;
	char filename[200];
#ifdef GROUPS
  char *Key = "Group";
#else
  char *Key = "Subhalo";
#endif
  enum particle_field {Gas, DarkMatter, Disk, Bulge, Stars, Bndry};
#ifdef GAS
  type_int type_particle[1] = {Gas};
  sprintf(filename,"%s_Gas",Key);
#elif  STARS
  type_int type_particle[1] = {Stars};
 	sprintf(filename,"%s_Stars",Key);
#elif STARS_GAS 
  type_int type_particle[2] = {Stars, Gas};
	sprintf(filename,"%s_Stars_Gas",Key);
#else
  type_int type_particle[1] = {DarkMatter};
	sprintf(filename,"%s_Dark_Matter",Key);
#endif
 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NumProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);

  if(MyProc == 0)
  {
    TIMER(start_tot);
    TIMER(start);
  }

  init_variables(argc,argv);

  MPI_Barrier(MPI_COMM_WORLD);

  leeheader(type_particle);

  MPI_Barrier(MPI_COMM_WORLD);

  read_gadget(type_particle, MyProc, NumProc);
  
  MPI_Barrier(MPI_COMM_WORLD);

  if(MyProc == 0)
  {
    TIMER(end);
    fprintf(stdout,"Total read time %f\n",end-start);
    fflush(stdout);
    TIMER(start);
  }

  correlacion(MyProc, NumProc, filename);
  free(tracers);

  MPI_Barrier(MPI_COMM_WORLD);

  if(MyProc == 0)
  {
    TIMER(end);
    TIMER(end_tot);
    fprintf(stdout,"Total Correlation time %f\n",end-start);
    fprintf(stdout,"Total time %f\n",end_tot-start_tot);
    fflush(stdout);
  }

  MPI_Finalize();

  return(EXIT_SUCCESS);
}
