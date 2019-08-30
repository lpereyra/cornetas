#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
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

static void auto_cross(const type_int NTHREADS, const type_int Ncentros, \
                       type_real * const centr, type_real * const tracer, type_real * const vdir, \
                       long unsigned *** restrict count_r, long unsigned * restrict Nc)
{
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
	type_int i, j, ivir, Tid;
	type_real delta[3];
	type_real r, direcc;
 	int i_r;
 
  #pragma omp parallel for schedule(static) num_threads(NTHREADS) \
  default(none) private(i, j, Tid, ix, iy, iz, ibox, \
  xe_sup, xe_inf, ye_inf, ye_sup, ze_sup, ze_inf, \
  delta,r,i_r,ivir,direcc) shared(grid, count_r, Nc, cp, stdout) 
  for(i=0; i<Ncentros; i++)
  {

		if(i%10000==0) 
    {
      fprintf(stdout,"%u/%u - %.2f\n",i,Ncentros,(float)i/(float)Ncentros);
      fflush(stdout);
		}

    Tid = omp_get_thread_num();
    Nc[Tid] += 1;

    ix = (long)((centr[Ndim*i])*fac);
    iy = (long)((centr[Ndim*i+1])*fac);
    iz = (long)((centr[Ndim*i+2])*fac);

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
            delta[0] = tracer[Ndim*j]   - centr[Ndim*i]  ; 
          	delta[1] = tracer[Ndim*j+1] - centr[Ndim*i+1];
          	delta[2] = tracer[Ndim*j+2] - centr[Ndim*i+2];
          
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
            
            count_r[Tid][Nvdir][i_r] += 1; // Isotropica

            for(ivir=0;ivir<Nvdir;ivir++)
            {
              //direcc  = fabs(vdir[ivir*Ncentros+Ndim*i]*delta[0] + vdir[ivir*Ncentros+Ndim*i+1]*delta[1] + vdir[ivir*Ncentros+Ndim*i+2]*delta[2]);
              direcc  = fabs(vdir[Ndim*(Nvdir*i+ivir)]*delta[0] + vdir[Ndim*(Nvdir*i+ivir)+1]*delta[1] + vdir[Ndim*(Nvdir*i+ivir)+2]*delta[2]);

              if(direcc > COS_ANGLE_CUT_PAR)      // PARALELA
                count_r[Tid][ivir][i_r] += 1;
              else if(direcc < COS_ANGLE_CUT_PER) // PERPENDICULAR
                count_r[Tid][ivir][i_r + NBINS] += 1;
            }

		  	  }
		    }		 
  		}	
    }
	}

  return;
}

#ifdef LANDY
  static void estimator(const type_int ivir, FILE *f_Xi_r_par, FILE *f_Xi_r_per, double *DD_r,  double *DR_r, double *RR_r,  double *RD_r)
#else
  static void estimator(const type_int ivir, FILE *f_Xi_r_par, FILE *f_Xi_r_per, double *DD_r,  double *DR_r)
#endif
{
  int i;
  type_real *Xi_r;

  if(ivir != Nvdir)
    Xi_r = (type_real *) malloc(2*NBINS*sizeof(type_real));
  else
    Xi_r = (type_real *) malloc(NBINS*sizeof(type_real));

  for(i=0;i<(ivir!=Nvdir ? (2*NBINS) : NBINS);i++)
  {
    #ifdef LANDY
      Xi_r[i] = RR_r[i]>0.0f ? (( DD_r[i] - DR_r[i] - RD_r[i] +   RR_r[i]) ) / ( RR_r[i]) : 0.0f;	//Landy & Szalay (1993)
    #else
      Xi_r[i] = DR_r[i]>0.0f ? (DD_r[i] -  DR_r[i]) / DR_r[i] : 0.0f;
    #endif
  }

  for(i=0;i<NBINS;i++)
  {
    fprintf(f_Xi_r_par,"%lf ",Xi_r[i]);
  }
  fprintf(f_Xi_r_par,"\n");
  
  if(ivir != Nvdir)
  {
    for(i=NBINS;i<(2*NBINS);i++)
    {
      fprintf(f_Xi_r_per,"%lf ",Xi_r[i]);
    }
    fprintf(f_Xi_r_per,"\n");
  }
  
  free(Xi_r);    

  return;

}

#ifdef LANDY
static void asignacion(type_real *rd_cent, type_real *vdir_rd_cent)
{
  long i;
  srand(80); // guarda la semilla
  type_real alfa, delta;

  for(i=0;i<fc.nrand_cent;i++)
  {
    rd_cent[Ndim*i]   = drand48()*cp.lbox;
    rd_cent[Ndim*i+1] = drand48()*cp.lbox;
    rd_cent[Ndim*i+2] = drand48()*cp.lbox;

    alfa  = drand48()*2.0f*M_PI;
    delta = 2.0*drand48()-1.0f;
    delta = asin(delta);

    vdir_rd_cent[Ndim*i]   = cos(delta)*cos(alfa);
    vdir_rd_cent[Ndim*i+1] = cos(delta)*sin(alfa);
    vdir_rd_cent[Ndim*i+2] = sin(delta);    
  }

  return;
}
#endif

extern void correlacion(void)
{
  type_int i, j, k, ivir, isample;
  long unsigned Nc_DD, Nc_DR;
#ifdef LANDY
  long unsigned Nc_RD, Nc_RR;
#endif
  const type_int NTHREADS = omp_get_max_threads();
  //const type_int NTHREADS = 26;
  const type_real cte_vol_par = (4./3.)*M_PI*(1.0 - COS_ANGLE_CUT_PAR); // EN EL CASO DE 60 grados
  const type_real cte_vol_per = (4./3.)*M_PI*COS_ANGLE_CUT_PER;         // AMBOS VOL COINCIDEN
  type_real *bins;
  long unsigned   *aux_Nc_DD, *aux_Nc_DR;
  long unsigned ***aux_DD_r;
  double *DD_r, *DR_r;
#ifdef LANDY
  long unsigned   *aux_Nc_RD, *aux_Nc_RR;
  long unsigned ***aux_RD_r;
  double *RD_r, *RR_r;
#endif
  FILE *fout;
  FILE *f_Xi_r_par;
  FILE *f_Xi_r_per;
  FILE *f_cont_DD_r;
  FILE *f_cont_DR_r;
#ifdef LANDY
  FILE *f_cont_RR_r;
  FILE *f_cont_RD_r;
#endif
  struct stat st = {0};
  char   folder[200];
  double start,end;
#ifdef SAMPLE_CONTROL
  char   texto[Nvdir+1][50] = {"momento_angular_","semieje_mayor_","semieje_intermedio_","semieje_menor_","isotropica_"};
#else
  char   texto[Nvdir+1][50] = {"momento_angular_","semieje_mayor_","semieje_intermedio_","semieje_menor_","filament_","isotropica_",};
#endif
#ifdef HALO_PARTICULA
  char   buff[200];
#endif
#ifdef LANDY
  type_real *rd_cent;
  type_real *vdir_rd_cent;
#endif

  for(isample=0;isample<NSAMPLES;isample++)
  {
    BLUE("*****************************\n");
    BLUE("*****************************\n");
    sprintf(message,"********** SAMPLE %u *********\n",isample);
    RED(message);
    BLUE("*****************************\n");
    BLUE("*****************************\n");
    fflush(stdout);

    read_grup_sample(fof, isample);

  #ifdef SAMPLE_CONTROL
     sprintf(folder,"out_cross_control_%s_ang_%.2f_Nrand_%.4d/",snap.sample_info,(type_real)ANGLE_CUT,NRAND);
  #else
     sprintf(folder,"out_cross_filament_%s_ang_%.2f_Nrand_%.4d/",snap.sample_info,(type_real)ANGLE_CUT,NRAND);
  #endif

  #ifdef HALO_PARTICULA
    sprintf(buff,"halo_particula_%.2f_percent_",100.0*FRACC);
    strcpy(folder,concat(buff,folder));
  #endif

  if(stat(folder, &st) == -1) 
  {
     mkdir(folder, 0700);
  }

  fprintf(stdout,"%s\n",folder);

  BLUE("*****************************\n");
  RED("********* ALLOC MEM *********\n");
  RED("******** INICIA CONT ********\n");
  BLUE("*****************************\n");

  fout   = fopen(concat(folder,"bins_cross.out"),"w");
  bins   = (type_real *) malloc((NBINS+1)*sizeof(type_real));
  init_bins(NBINS+1,bins);
  for(i=0;i<(NBINS+1);i++) 
    fprintf(fout,"%.16f ",bins[i]);
  fprintf(fout,"\n");
  fclose(fout);
  
  DD_r = (double *) calloc(2*NBINS,sizeof(double));
  DR_r = (double *) calloc(2*NBINS,sizeof(double));

  aux_Nc_DD = (long unsigned   *) calloc(NTHREADS,sizeof(long unsigned  ));
  aux_Nc_DR = (long unsigned   *) calloc(NTHREADS,sizeof(long unsigned  ));
  aux_DD_r  = (long unsigned ***) calloc(NTHREADS,sizeof(long unsigned **));

#ifdef LANDY
  RD_r = (double *) calloc(2*NBINS,sizeof(double));
  RR_r = (double *) calloc(2*NBINS,sizeof(double));

  aux_Nc_RD = (long unsigned   *) calloc(NTHREADS,sizeof(long unsigned  ));
  aux_Nc_RR = (long unsigned   *) calloc(NTHREADS,sizeof(long unsigned  ));
  aux_RD_r  = (long unsigned ***) calloc(NTHREADS,sizeof(long unsigned **));
#endif

  for(i=0;i<NTHREADS;i++)
  {
    aux_DD_r[i] = (long unsigned **) calloc((Nvdir+1),sizeof(long unsigned *));
#ifdef LANDY    
    aux_RD_r[i] = (long unsigned **) calloc((Nvdir+1),sizeof(long unsigned *));
#endif
    for(ivir=0;ivir<=Nvdir;ivir++)
    {
      aux_DD_r[i][ivir] = (long unsigned *) calloc((ivir!=Nvdir ? (2*NBINS) : NBINS),sizeof(long unsigned));
  #ifdef LANDY   
      aux_RD_r[i][ivir] = (long unsigned *) calloc((ivir!=Nvdir ? (2*NBINS) : NBINS),sizeof(long unsigned));
  #endif
    }
  }

  BLUE("*****************************\n");
  BLUE("*****************************\n");
  fflush(stdout);

  RED("********* ALLOC MEM *********\n");
  RED("******* INICIA ASSIG ********\n");
  BLUE("*****************************\n");

  fc.ndata_cent = (long unsigned)cp.ngrup_sample;
  fc.nrand_cent = (long unsigned)cp.ngrup_sample*NRAND;

#ifdef HALO_PARTICULA
  fc.ndata_trac  = (long unsigned)cp.npart;
  fc.nrand_trac  = (long unsigned)cp.npart*NRAND;
#else
  fc.ndata_trac  = (long unsigned)cp.ngrup;
  fc.nrand_trac  = (long unsigned)cp.ngrup*NRAND;
#endif

#ifdef LANDY
  rd_cent = (type_real *) malloc(Ndim*fc.nrand_cent*sizeof(type_real));
  vdir_rd_cent = (type_real *) malloc(Ndim*fc.nrand_cent*sizeof(type_real));
  asignacion(rd_cent,vdir_rd_cent);
#endif

  fout   = fopen(concat(folder,"ntotal.out"),"w");
  fprintf(fout,"#ncentr nrand_cent ntrac nrand_trac\n");
  fprintf(fout,"%lu %lu %lu %lu\n",fc.ndata_cent,fc.nrand_cent,fc.ndata_trac,fc.nrand_trac);
  fclose(fout);

  BLUE("*****************************\n");
  BLUE("*****************************\n");
  fflush(stdout);
 
  RED("***** INICIA CORRELACION ****\n");
  BLUE("*****************************\n");
  RED("****** INICIA DAT-DAT *******\n");
  BLUE("*****************************\n");

  if(isample==0) grid_init(fc.ndata_trac, tracers, END_BINS);

  TIMER(start);
  auto_cross(NTHREADS,fc.ndata_cent,centros,tracers,vdir_centros,aux_DD_r,aux_Nc_DD);
  TIMER(end);

  for(ivir=0;ivir<=Nvdir;ivir++)
  {
    f_cont_DD_r = fopen(concat(folder,concat(texto[ivir],"threads_contadores_DD_r.out")),"w");
    for(i=0;i<NTHREADS;i++)
    {
      fprintf(f_cont_DD_r,"%lu ",aux_Nc_DD[i]);
      for(j=0;j<(ivir!=Nvdir ? (2*NBINS) : NBINS);j++)
        fprintf(f_cont_DD_r,"%lf ",(double)aux_DD_r[i][ivir][j]);
      fprintf(f_cont_DD_r,"\n");
    }
    fclose(f_cont_DD_r);
  }

  fprintf(stdout,"Total time DAT-DAT %f\n",end-start);
  fflush(stdout);

  BLUE("*****************************\n");
  RED("*****************************\n");

  for(ivir=0;ivir<=Nvdir;ivir++)
  {
    f_cont_DR_r = fopen(concat(folder,concat(texto[ivir],"threads_contadores_DR_r.out")),"w");
    for(i=0;i<NTHREADS;i++)
    {
      aux_Nc_DR[i] = aux_Nc_DD[i];
      
      fprintf(f_cont_DR_r,"%lu ",aux_Nc_DR[i]);

      for(j=0;j<NBINS;j++)
      {
        fprintf(f_cont_DR_r,"%lf ",(double)aux_Nc_DR[i]*\
        (ivir != Nvdir ?  cte_vol_par : (4./3.)*M_PI)*\
        (pow(bins[j+1],3)-pow(bins[j],3))*\
        ((double)fc.nrand_trac/pow(cp.lbox,3)));
      }

      for(j=0;j<(ivir!=Nvdir ? NBINS : 0);j++)
      {
        fprintf(f_cont_DR_r,"%lf ",(double)aux_Nc_DR[i]*\
        cte_vol_per*(pow(bins[j+1],3)-pow(bins[j],3))*\
        ((double)fc.nrand_trac/pow(cp.lbox,3)));
      }
      
      fprintf(f_cont_DR_r,"\n");
    }
    fclose(f_cont_DR_r); 
  }

#ifdef LANDY

  BLUE("*****************************\n");
  RED("****** INICIA RAN-DAT *******\n");
  BLUE("*****************************\n");

  TIMER(start);
  auto_cross(NTHREADS,fc.nrand_cent,rd_cent,tracers,vdir_rd_cent,aux_RD_r,aux_Nc_RD);
  TIMER(end);

  for(ivir=0;ivir<=Nvdir;ivir++)
  {
    f_cont_RD_r = fopen(concat(folder,concat(texto[ivir],"threads_contadores_RD_r.out")),"w");
    for(i=0;i<NTHREADS;i++)
    {
      fprintf(f_cont_RD_r,"%lu ",aux_Nc_RD[i]);
      for(j=0;j<(ivir!=Nvdir ? (2*NBINS) : NBINS);j++)
        fprintf(f_cont_RD_r,"%lf ",(double)aux_RD_r[i][ivir][j]);
      fprintf(f_cont_RD_r,"\n");
    }
    fclose(f_cont_RD_r);
  }

  fprintf(stdout,"Total time RAN-DAT %f\n",end-start);
  fflush(stdout);

  BLUE("*****************************\n");
  RED("*****************************\n");

  for(ivir=0;ivir<=Nvdir;ivir++)
  {
    f_cont_RR_r = fopen(concat(folder,concat(texto[ivir],"threads_contadores_RR_r.out")),"w");

    for(i=0;i<NTHREADS;i++)
    {
      aux_Nc_RR[i] = aux_Nc_RD[i];
 
      fprintf(f_cont_RR_r,"%lu ",aux_Nc_RR[i]);

      for(j=0;j<NBINS;j++)
      {
        fprintf(f_cont_RR_r,"%lf ",(double)aux_Nc_RR[i]*\
        (ivir != Nvdir ?  cte_vol_par : (4./3.)*M_PI)*\
        (pow(bins[j+1],3)-pow(bins[j],3))*\
        ((double)fc.nrand_trac/pow(cp.lbox,3)));
      }

      for(j=0;j<(ivir!=Nvdir ? NBINS : 0);j++)
      { 
        fprintf(f_cont_RR_r,"%lf ",(double)aux_Nc_RR[i]*\
        cte_vol_per*(pow(bins[j+1],3)-pow(bins[j],3))*\
        ((double)fc.nrand_trac/pow(cp.lbox,3)));
      }

      fprintf(f_cont_RR_r,"\n");
    }
    fclose(f_cont_RR_r); 
  }

#endif

  for(ivir=0;ivir<=Nvdir;ivir++)
  {
    if(ivir!=Nvdir)
    {
      f_Xi_r_par = fopen(concat(folder,concat(texto[ivir],"Xi_r_cross_par.out")),"w");
      f_Xi_r_per = fopen(concat(folder,concat(texto[ivir],"Xi_r_cross_per.out")),"w");
    }else{
      f_Xi_r_par = fopen(concat(folder,concat(texto[ivir],"Xi_r_cross.out")),"w");     
      f_Xi_r_per = NULL;
    }
  
    for(k=0;k<(NTHREADS+1);k++)
    {    
  
  #ifdef LANDY
      Nc_DD = Nc_DR = Nc_RD = Nc_RR = 0;
  #else
      Nc_DD = Nc_DR = 0;
  #endif
      for(i=0;i<NTHREADS;i++)
      {
        if(i==k) continue;
  
        Nc_DD += aux_Nc_DD[i];
        Nc_DR += aux_Nc_DR[i];
  #ifdef LANDY
        Nc_RD += aux_Nc_RD[i];
        Nc_RR += aux_Nc_RR[i];
  #endif
      }
  
      for(j=0;j<(2*NBINS);j++)
  #ifdef LANDY
        DD_r[j] = DR_r[j] = RR_r[j] = RD_r[j] = 0.0f;
  #else
        DD_r[j] = DR_r[j] = 0.0f;
  #endif
  
      for(i=0;i<NTHREADS;i++)
      {
        if(i==k) continue;
  
        for(j=0;j<NBINS;j++)
        {
          DD_r[j] += ((double)aux_DD_r[i][ivir][j] / (double)Nc_DD / (double)fc.ndata_trac);
          DR_r[j] += (double)aux_Nc_DR[i]*(ivir != Nvdir ?  cte_vol_par : (4./3.)*M_PI)*(pow(bins[j+1],3)-pow(bins[j],3)) \
                     * ((double)fc.nrand_trac/pow(cp.lbox,3)) \
                     / (double)Nc_DR / (double)fc.nrand_trac;
      #ifdef LANDY
          RD_r[j] += ((double)aux_RD_r[i][ivir][j] / (double)Nc_RD / (double)fc.ndata_trac);
          RR_r[j] += (double)aux_Nc_RR[i]*(ivir != Nvdir ?  cte_vol_par : (4./3.)*M_PI)*(pow(bins[j+1],3)-pow(bins[j],3)) \
                     * ((double)fc.nrand_trac/pow(cp.lbox,3)) \
                     / (double)Nc_RR / (double)fc.nrand_trac;
      #endif
        }
  
        for(j=NBINS;j<(ivir!=Nvdir ? (2*NBINS) : NBINS);j++)
        {
          DD_r[j] += ((double)aux_DD_r[i][ivir][j] / (double)Nc_DD / (double)fc.ndata_trac);
          DR_r[j] += (double)aux_Nc_DR[i]*cte_vol_per*(pow(bins[j-NBINS+1],3)-pow(bins[j-NBINS],3)) \
                     * ((double)fc.nrand_trac/pow(cp.lbox,3)) \
                     / (double)Nc_DR / (double)fc.nrand_trac;
      #ifdef LANDY
          RD_r[j] += ((double)aux_RD_r[i][ivir][j] / (double)Nc_RD / (double)fc.ndata_trac);
          RR_r[j] += (double)aux_Nc_RR[i]*cte_vol_per*(pow(bins[j-NBINS+1],3)-pow(bins[j-NBINS],3))* \
                     * ((double)fc.nrand_trac/pow(cp.lbox,3)) \
                     / (double)Nc_RR / (double)fc.nrand_trac;
      #endif
        }
      }
  
      #ifdef LANDY
        estimator(ivir, f_Xi_r_par, f_Xi_r_per, DD_r, DR_r, RR_r, RD_r);
      #else 
        estimator(ivir, f_Xi_r_par, f_Xi_r_per, DD_r, DR_r);
      #endif
    }
  
    fclose(f_Xi_r_par);
    if(ivir != Nvdir)
      fclose(f_Xi_r_per);
  }

  BLUE("*****************************\n");

  for(i=0;i<NTHREADS;i++)
  {
    for(ivir=0;ivir<=Nvdir;ivir++)
    {
       free(aux_DD_r[i][ivir]);
      #ifdef LANDY
        free(aux_RD_r[i][ivir]);
      #endif
    }

    free(aux_DD_r[i]);
    #ifdef LANDY
      free(aux_RD_r[i]);
    #endif
  }

  free(DD_r);
  free(DR_r);
  free(aux_DD_r);
  free(aux_Nc_DD);
  free(aux_Nc_DR);
#ifdef LANDY
  free(RD_r);
  free(RR_r);
  free(aux_RD_r);
  free(aux_Nc_RD);
  free(aux_Nc_RR);
#endif
  free(bins);

  free(centros);
  free(vdir_centros);
#ifdef LANDY
  free(rd_cent);
  free(vdir_rd_cent);
#endif
  if(isample==NSAMPLES-1)
    grid_free();
  
  } // isample

  return;

}
