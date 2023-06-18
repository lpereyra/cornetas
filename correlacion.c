#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "variables.h"
#include "cosmoparam.h"
#include "correlacion.h"
#include "leefile.h"
#include "colores.h"
#include "timer.h"
#include "grid.h"

#define IDX(j,i,n) j*n+i
static const type_real COS_ANGLE_CUT_PAR = cos(ANGLE_CUT*M_PI/180.0f);
static const type_real COS_ANGLE_CUT_PER = 1.0f - cos(ANGLE_CUT*M_PI/180.0f);

struct gridst grid;

static char *concat(const char *s1, const char *s2)
{
  char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
  // in real code you would check for errors in malloc here
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}

static void init_bins(const type_int nbins, type_real *bins)
{
  type_int i;
	type_real	power;
#ifdef LOG
	type_real	base_min = log10(START_BINS);
	type_real	step     = (log10(END_BINS)-log10(START_BINS))/(type_real)(nbins-1);
#else
	type_real	base_min = START_BINS;
	type_real	step     = (END_BINS-START_BINS)/(type_real)(nbins-1);
#endif

	for (i=0; i<nbins; i++)
  {
		power = base_min+(type_real)i*step;
    #ifdef LOG
		  bins[i] = pow(10.0,power);
    #else
		  bins[i] = power;
    #endif
  }

  return;
}

inline static long igrid(const long ix, const long iy, const long iz, const long ngrid)
{
  return (ix * ngrid + iy) * ngrid + iz;
}

static void auto_cross(const type_int Ncentros, \
                       type_real * const centr, type_real * const tracer, type_real * const vdir, \
                       long unsigned *** restrict count_r)
{
  const type_int NTHREADS = omp_get_max_threads();
#ifdef LOG
  const type_real start_bin = log10(START_BINS);
  const type_real end_bin   = log10(END_BINS);
#else
  const type_real start_bin = START_BINS;
  const type_real end_bin   = END_BINS;
#endif
  const type_real step = (end_bin-start_bin)/(type_real)(NBINS);
  const type_real fac = (type_real)grid.ngrid/cp.lbox;
  long xe_sup, xe_inf, ye_sup, ye_inf, ze_sup, ze_inf;
  long ix, iy, iz, ibox;
  type_int i, j, ivir;
  type_real delta[3];
  type_real r, direcc;
  int i_r;

#ifdef __GNUC__
  // If  gcc_version >= 10.0
  #  if __GNUC_PREREQ(10,0)
    #pragma omp parallel for schedule(static) num_threads(NTHREADS) \
    default(none) private(i, j, ix, iy, iz, ibox, \
    xe_sup, xe_inf, ye_inf, ye_sup, ze_sup, ze_inf, \
    delta,r,i_r,ivir,direcc) shared(grid, count_r, cp, fac, \
    start_bin, end_bin, step, Ncentros, centr, tracer, vdir, stdout)
  #else
    #pragma omp parallel for schedule(static) num_threads(NTHREADS) \
    default(none) private(i, j, ix, iy, iz, ibox, \
    xe_sup, xe_inf, ye_inf, ye_sup, ze_sup, ze_inf, \
    delta,r,i_r,ivir,direcc) shared(grid, count_r, cp, stdout) 
  #endif
#endif
  for(i=0; i<Ncentros; i++)
  {

		if(i%10000==0) 
    {
      fprintf(stdout,"%u/%u - %.2f\n",i,Ncentros,(float)i/(float)Ncentros);
      fflush(stdout);
		}

    ix = (long)((centr[Ndim*i])*fac);
    iy = (long)((centr[Ndim*i+1])*fac);
    iz = (long)((centr[Ndim*i+2])*fac);

#ifdef SKIP
    if(ix < 0 && iy < 0 && iz < 0)
      continue;
#endif

    xe_inf = ix-1;
    xe_sup = ix+1;
              
    ye_inf = iy-1;
    ye_sup = iy+1;
              
    ze_inf = iz-1;
    ze_sup = iz+1;
                                                                  
    #ifndef PERIODIC
		if(xe_inf < 0) xe_inf = 0; 
		if(ye_inf < 0) ye_inf = 0; 
		if(ze_inf < 0) ze_inf = 0; 
    if(xe_inf >= grid.ngrid) xe_inf = (long)grid.ngrid-1;
    if(ye_inf >= grid.ngrid) ye_inf = (long)grid.ngrid-1;
    if(ze_inf >= grid.ngrid) ze_inf = (long)grid.ngrid-1;

    if(xe_sup < 0) xe_sup = 0; 
		if(ye_sup < 0) ye_sup = 0; 
		if(ze_sup < 0) ze_sup = 0; 
    if(xe_sup >= grid.ngrid) xe_sup = (long)grid.ngrid-1;
    if(ye_sup >= grid.ngrid) ye_sup = (long)grid.ngrid-1;
    if(ze_sup >= grid.ngrid) ze_sup = (long)grid.ngrid-1;
    #endif

    for(iy=ye_inf; iy<=ye_sup; iy++)
    {
      for(ix=xe_inf; ix<=xe_sup; ix++)
      {
		    for(iz=ze_inf; iz<=ze_sup; iz++)
        {
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

          for(j=grid.start[ibox];j<( ibox == grid.nalloc-1 ? grid.nobj : grid.start[ibox+1]);j++)
          {
            delta[0] = tracer[(long unsigned)Ndim*j]   - centr[Ndim*i]  ; 
          	delta[1] = tracer[(long unsigned)Ndim*j+1] - centr[Ndim*i+1];
          	delta[2] = tracer[(long unsigned)Ndim*j+2] - centr[Ndim*i+2];
          
            #ifdef PERIODIC
            delta[0] = ( delta[0] >  cp.lbox*0.5f ) ? delta[0] - cp.lbox : delta[0] ;
            delta[1] = ( delta[1] >  cp.lbox*0.5f ) ? delta[1] - cp.lbox : delta[1] ;
            delta[2] = ( delta[2] >  cp.lbox*0.5f ) ? delta[2] - cp.lbox : delta[2] ;
            delta[0] = ( delta[0] < -cp.lbox*0.5f ) ? delta[0] + cp.lbox : delta[0] ;
            delta[1] = ( delta[1] < -cp.lbox*0.5f ) ? delta[1] + cp.lbox : delta[1] ;
            delta[2] = ( delta[2] < -cp.lbox*0.5f ) ? delta[2] + cp.lbox : delta[2] ;
            #endif
          
            r = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
          
          	if((r   <= START_BINS) || (r   >= END_BINS)) continue;
           
            delta[0] /= r;  
            delta[1] /= r;  
            delta[2] /= r;

            #ifdef LOG
            r = log10(r);
            #endif
          
          	i_r = (int)((r - start_bin)/step);
          
            if(i_r == NBINS) i_r = NBINS-1;
            if(i_r == -1)    i_r = 0;                  
            
            count_r[i][Nvdir][i_r] += 1; // Isotropica

            for(ivir=0;ivir<Nvdir;ivir++)
            {
              //direcc  = fabs(vdir[ivir*Ncentros+Ndim*i]*delta[0] + vdir[ivir*Ncentros+Ndim*i+1]*delta[1] + vdir[ivir*Ncentros+Ndim*i+2]*delta[2]);
              direcc  = fabs(vdir[Ndim*(Nvdir*i+ivir)]*delta[0] + vdir[Ndim*(Nvdir*i+ivir)+1]*delta[1] + vdir[Ndim*(Nvdir*i+ivir)+2]*delta[2]);

              if(direcc > COS_ANGLE_CUT_PAR)      // PARALELA
                count_r[i][ivir][i_r] += 1;
              else if(direcc < COS_ANGLE_CUT_PER) // PERPENDICULAR
                count_r[i][ivir][i_r + NBINS] += 1;
            }

		  	  }
		    }		 
  		}	
    }
	}

  return;
}

extern void correlacion(const int MyProc, const int NumProc, const char *filetype)
{
  type_int i, j, ivir;
  type_real *bins;
  FILE *fout;
  FILE *f_cont_DD_r;
  long unsigned ***aux_DD_r;
  long unsigned ***aux_DD_r_global;
  struct stat st = {0};
  char   folder[200];
  double start,end;
  char   texto[Nvdir+1][50] = {"momento_angular_","semieje_mayor_","semieje_intermedio_","semieje_menor_","isotropica_",};
  char   buff[200];
  char   filename[200];
  struct fcorrST fc;
  long unsigned aux_global;

  read_grup_sample(MyProc, NumProc, filetype);

  sprintf(folder,"total_out_cross_ang_%.2f_Nrand_%.4d/",(type_real)ANGLE_CUT,NRAND);

  #ifdef SKIP
    sprintf(buff,"%s_%.2f_%.2f_halo_particula_",filetype,MCUT_MIN,MCUT_MAX);
  #else
    sprintf(buff,"%s_halo_particula_",filetype);
  #endif
  strcpy(folder,concat(buff,folder));

  if(MyProc == 0)
  {
    if(stat(folder, &st) == -1) 
    {
       mkdir(folder, 0700);
    }
    fprintf(stdout,"%s\n",folder);
  }

  bins   = (type_real *) malloc((NBINS+1)*sizeof(type_real));
  init_bins(NBINS+1,bins);

  if(MyProc == 0)
  {
    sprintf(filename,concat(folder,"bins_cross.out"));
    fout   = fopen(filename,"w");
    for(i=0;i<(NBINS+1);i++) 
      fprintf(fout,"%.16f ",bins[i]);
    fprintf(fout,"\n");
    fclose(fout);

    BLUE("*****************************\n");
    RED("********* ALLOC MEM *********\n");
    RED("******** INICIA CONT ********\n");
    BLUE("*****************************\n");
  }

  MPI_Barrier (MPI_COMM_WORLD);

  fc.ndata_cent = (long unsigned)cp.ngrup_sample;
  fc.nrand_cent = (long unsigned)cp.ngrup_sample*NRAND;
  fc.ndata_trac  = (long unsigned)cp.npart;
  fc.nrand_trac  = (long unsigned)cp.npart*NRAND;
 
  aux_DD_r         = (long unsigned ***) calloc(fc.ndata_cent,sizeof(long unsigned **));
  aux_DD_r_global  = (long unsigned ***) calloc(fc.ndata_cent,sizeof(long unsigned **));

  for(i=0;i<fc.ndata_cent;i++)
  {
    aux_DD_r[i]        = (long unsigned **) calloc((Nvdir+1),sizeof(long unsigned *));
    aux_DD_r_global[i] = (long unsigned **) calloc((Nvdir+1),sizeof(long unsigned *));
    for(ivir=0;ivir<=Nvdir;ivir++)
    {
      aux_DD_r[i][ivir]        = (long unsigned *) calloc((ivir!=Nvdir ? (2*NBINS) : NBINS),sizeof(long unsigned));
      aux_DD_r_global[i][ivir] = (long unsigned *) calloc((ivir!=Nvdir ? (2*NBINS) : NBINS),sizeof(long unsigned));
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if(MyProc == 0)
  {
    BLUE("*****************************\n");
    BLUE("*****************************\n");
    fflush(stdout);
 
    RED("***** INICIA CORRELACION ****\n");
    BLUE("*****************************\n");
    RED("****** INICIA DAT-DAT *******\n");
    BLUE("*****************************\n");
    fflush(stdout);
  }

  grid_init(fc.ndata_trac, tracers, END_BINS);
  TIMER(start);
  auto_cross(fc.ndata_cent,centros,tracers,vdir_centros,aux_DD_r);
  TIMER(end);

  MPI_Barrier (MPI_COMM_WORLD);

  if(MyProc == 0)
  {
    BLUE("*****************************\n");
    RED("***** REDUCE CONTADORES *****\n");
    fflush(stdout);
  }

  MPI_Reduce(&fc.ndata_trac, &aux_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  fc.ndata_trac = aux_global;
  MPI_Reduce(&fc.nrand_trac, &aux_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  fc.nrand_trac = aux_global;

  for(i=0;i<fc.ndata_cent;i++)
  {
    for(ivir=0;ivir<=Nvdir;ivir++)
    {
      j = ivir!=Nvdir ? (2*NBINS) : NBINS;
      MPI_Reduce(aux_DD_r[i][ivir], aux_DD_r_global[i][ivir], j, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if(MyProc == 0)
  {
    BLUE("*****************************\n");
    fprintf(stdout,"Total time DAT-DAT %f\n",end-start); fflush(stdout);
    BLUE("*****************************\n");
    
    RED("************ WRITE **********\n");
    BLUE("*****************************\n");
    fflush(stdout);
    
    TIMER(start);
    
    sprintf(filename,concat(folder,"ntotal.out"));
    fout   = fopen(filename,"w");
    fprintf(fout,"#ncentr nrand_cent ntrac nrand_trac\n");
    fprintf(fout,"%lu %lu %lu %lu\n",fc.ndata_cent,fc.nrand_cent,fc.ndata_trac,fc.nrand_trac);
    fclose(fout);

    for(ivir=0;ivir<=Nvdir;ivir++)
    {
      sprintf(filename,concat(folder,concat(texto[ivir],"threads_contadores_DD_r.out")));
      f_cont_DD_r = fopen(filename,"w");
      sprintf(filename,concat(folder,concat(texto[ivir],"threads_contadores_DD_r.bin")));
      fout        = fopen(filename,"wb");
      fwrite(&fc.ndata_cent,sizeof(long unsigned),1,fout);
      
      for(i=0;i<fc.ndata_cent;i++)
      {
        for(j=0;j<(ivir!=Nvdir ? (2*NBINS) : NBINS);j++)
        {
          fprintf(f_cont_DD_r,"%lu ",aux_DD_r_global[i][ivir][j]);
          fwrite(&aux_DD_r_global[i][ivir][j],sizeof(long unsigned),1,fout);
        }
        fprintf(f_cont_DD_r,"\n");
      }
      
      fclose(f_cont_DD_r);
      fclose(fout);
    }

    TIMER(end);
    
    fprintf(stdout,"Total Write %f\n",end-start); fflush(stdout);

    BLUE("*****************************\n");
    RED("*****************************\n");
    BLUE("*****************************\n");
  }

  MPI_Barrier (MPI_COMM_WORLD);

  for(i=0;i<fc.ndata_cent;i++)
  {
    for(ivir=0;ivir<=Nvdir;ivir++)
    {
      free(aux_DD_r[i][ivir]);
      free(aux_DD_r_global[i][ivir]);
    }
    free(aux_DD_r[i]);
    free(aux_DD_r_global[i]);
  }

  free(aux_DD_r);
  free(aux_DD_r_global);
  free(bins);

  free(centros);
  free(vdir_centros);
  grid_free();

  return;
}
