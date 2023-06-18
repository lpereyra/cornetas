#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.h"
#include "cosmoparam.h"
#include "correlacion.h"
#include "grid.h"

inline static long igrid(const long ix, const long iy, const long iz, const long ngrid)
{
  return (ix * ngrid + iy) * ngrid + iz;
}

static void grid_build(type_real *Pos)
{
  long long i, j;
  long long ix, iy, iz;
	long long ibox;
  double fac;
  type_int *Id;
  unsigned long tmp;
  type_real Pos_source[3], Pos_save[3];
	type_int idsource,dest,idsave;

  fac = (double)grid.ngrid/(double)cp.lbox ;
	fprintf(stdout,"Building Grid..... Ngrid = %ld\n",grid.ngrid);
  fflush(stdout);

  Id = (type_int *) malloc(grid.nobj*sizeof(type_int));
  assert(Id != NULL);

  for(i=0;i<(long long)grid.nobj;i++)
  {

    ix = (long long)((double)(Pos[((long unsigned)Ndim*i)])  *fac);
    iy = (long long)((double)(Pos[((long unsigned)Ndim*i)+1])*fac);
    iz = (long long)((double)(Pos[((long unsigned)Ndim*i)+2])*fac);

#ifdef PERIODIC
    ibox = igrid(( (ix >= (long)grid.ngrid) ? ix - (long)grid.ngrid : ( (ix<0) ? ix + (long)grid.ngrid : ix ) ),\
                 ( (iy >= (long)grid.ngrid) ? iy - (long)grid.ngrid : ( (iy<0) ? iy + (long)grid.ngrid : iy ) ),\
                 ( (iz >= (long)grid.ngrid) ? iz - (long)grid.ngrid : ( (iz<0) ? iz + (long)grid.ngrid : iz ) ),\
                 (long)grid.ngrid);
#else
    ibox = igrid(( (ix >= (long)grid.ngrid) ? (long)grid.ngrid : ( (ix<0) ? 0 : ix ) ),\
                 ( (iy >= (long)grid.ngrid) ? (long)grid.ngrid : ( (iy<0) ? 0 : iy ) ),\
                 ( (iz >= (long)grid.ngrid) ? (long)grid.ngrid : ( (iz<0) ? 0 : iz ) ),\
                 (long)grid.ngrid);
#endif

    assert(ibox>=0 && ibox<((long long)grid.ngrid*(long long)grid.ngrid*(long long)grid.ngrid));

    Id[i] = grid.start[ibox];
    grid.start[ibox] = i;
  }

	fprintf(stdout,"FINISH first part Building Grid\n");
  fflush(stdout);

  j = 0;
  for(ibox=0;ibox<grid.nalloc;ibox++)
  {      
    i = grid.start[ibox];
    grid.start[ibox] = j;
    while(i != grid.nobj)
    {
      tmp = Id[i];
      Id[i] = j;
      i = tmp;
      j++;	
    }
  }

  fprintf(stdout,"BEGING reorder part Building Grid\n");
  fflush(stdout);

  for(i=0;i<grid.nobj;i++)
  {
    if(Id[i] != i)
    {
      Pos_source[0] = Pos[(long unsigned)Ndim*i];
      Pos_source[1] = Pos[(long unsigned)Ndim*i+1];
      Pos_source[2] = Pos[(long unsigned)Ndim*i+2];
    	idsource = Id[i];
      dest = Id[i];

      while(1)
      {
        Pos_save[0] = Pos[(long unsigned)Ndim*dest];
        Pos_save[1] = Pos[(long unsigned)Ndim*dest+1];
        Pos_save[2] = Pos[(long unsigned)Ndim*dest+2];
	      idsave = Id[dest];

        Pos[(long unsigned)Ndim*dest]   = Pos_source[0];
        Pos[(long unsigned)Ndim*dest+1] = Pos_source[1];
        Pos[(long unsigned)Ndim*dest+2] = Pos_source[2];
        Id[dest] = idsource;

        if(dest == i)  break;

        Pos_source[0] = Pos_save[0];
        Pos_source[1] = Pos_save[1];
        Pos_source[2] = Pos_save[2];
        idsource = idsave;
  	    dest = idsource;
      } // cierra el while
    }  // cierra el if
  } // cierra el for

  free(Id);

  fflush(stdout);

  return;
}

extern void grid_init(type_int nobj, type_real *Pos, const type_int cell_size)
{
  long unsigned i;

  grid.nobj = nobj;
  grid.ngrid = (long)(cp.lbox/cell_size);

  if(grid.ngrid > NGRIDMAX)
  {
    grid.ngrid = NGRIDMAX;
    fprintf(stdout,"Using NGRIDMAX = %ld\n",(long)NGRIDMAX);
  }else{
    fprintf(stdout,"Using NGRID    = %ld\n",grid.ngrid);
  }

  grid.nalloc = grid.ngrid*grid.ngrid*grid.ngrid;

  fprintf(stdout,"allocating %.5f gb\n",
                     (double)((grid.nalloc)*sizeof(type_int))/1024.0/1024.0/1024.0);

  grid.start	= (type_int *) malloc(grid.nalloc*sizeof(type_int));
  assert(grid.start != NULL);

  // Inicializa el arreglo
  for(i=0;i<grid.nalloc;i++)
    grid.start[i] = grid.nobj;

  grid_build(Pos);

  fflush(stdout);

  return;
}

extern void grid_free(void)
{
  if(grid.start!=NULL) free(grid.start);
  return;
}
