#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef NRAND
  #define NRAND 100
#endif

#ifndef NBINS
  #define NBINS         30     // NUMERO DE BINES
#endif

#ifndef START_BINS
  #define START_BINS  0.10*POSFACTOR
#endif

#ifndef END_BINS
  #define END_BINS    50.0*POSFACTOR
#endif

#ifndef ANGLE_CUT
  #define ANGLE_CUT 60 //M_PI
#endif

#ifndef TYPE_FLAG
  #define TYPE_FLAG 2           // Tipo dos
#endif

#ifdef CUT_IN_LEN

  #ifndef LEN_MIN
    #define LEN_MIN 10000         // en Kpc
  #endif

  #ifndef LEN_MAX
    #define LEN_MAX 100000        // en Kpc
  #endif

#endif

#ifdef PERCENT
  #ifndef FRACC 
    #define FRACC   0.1           //  EL PORCENTAJE QUE GUARDO
  #endif
#endif

#define NSAMPLES     6    /* Number of samples        */
#define N_part_types 6    /* Number of particle types */
#define Ndim         3    /* Number of dimension      */
#ifdef SAMPLE_CONTROL  
  #define Nvdir        4    /* Number of vdir           */
#else
  #define Nvdir        5    /* Number of vdir           */
#endif

#define START_BINS  0.10*POSFACTOR
#define END_BINS    50.0*POSFACTOR

/* Precision del codigo (reales) */
#ifdef PRECDOUBLE
typedef double type_real;
#else
typedef float type_real;
#endif

/* Precision del codigo (enteros) */
#ifdef LONGIDS
typedef unsigned long long type_int;
#else
typedef unsigned int type_int;
#endif

//struct grup_data
//{
//  type_int   save;
//  type_int   id;
//  type_real  Pos[3];
//  type_int   NumPart;
//  type_real  Mass;
//  type_real  aa;
//  type_real  bb;
//  type_real  cc;
//  type_real  ParP;
//  type_real  vcm[3];
//  type_real  L[3];
//  type_real  evec[9];
//  type_real  vdir[3];
//};

//#ifdef HALO_PARTICULA
//  struct particle_data 
//  {
//    //type_real      *x;
//    //type_real      *y;
//    //type_real      *z;
//    //#ifdef STORE_VELOCITIES
//    //type_real      *vx;
//    //type_real      *vy;
//    //type_real      *vz;
//    //#endif
//   
//    type_real Pos[3];
//    #ifdef STORE_VELOCITIES
//    type_real Vel[3];
//    #endif
//    #ifdef STORE_IDS
//    type_int       *id;
//    #endif
//  };
//#endif

extern type_int  nfrac;
extern type_real *fof;
#ifdef MCRITIC
extern type_real m_critica;
#endif
extern type_int  NNN;

extern type_real *tracers;
extern type_real *centros;
extern type_real *vdir_centros;

extern void init_variables(int argc, char **argv);

#endif
