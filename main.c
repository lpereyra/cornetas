#include <stdio.h>
#include <stdlib.h>

#include "variables.h"
#include "cosmoparam.h"
#include "correlacion.h"
#include "leefile.h"
#include "timer.h"

int main(int argc, char **argv)
{
  double start,end;
  double start_tot,end_tot;
  TIMER(start_tot);

  init_variables(argc,argv);
  NNN  = atoi(argv[2]);

  TIMER(start);

  leeheader();
  //read_grup_sample(fof);
#ifdef HALO_PARTICULA
  read_gadget();
#else
  read_grup_fof(fof);
#endif
  TIMER(end);
  fprintf(stdout,"Total read time %f\n",end-start);
  fflush(stdout);
  
  TIMER(start);
  correlacion();
  TIMER(end);
  fprintf(stdout,"Total correlation time %f\n",end-start);
  fflush(stdout);
  
  free(tracers);
  //free(centros);
  //free(vdir_centros);

  TIMER(end_tot);
  fprintf(stdout,"Total time %f\n",end_tot-start_tot);
  fflush(stdout);

  return(EXIT_SUCCESS);
}
