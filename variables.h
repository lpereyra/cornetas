#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#if defined(STARS_GAS)
  #define NTYPES 2
#else
  #define NTYPES 1
#endif

#ifndef NRAND
  #define NRAND 100
#endif

#ifndef NBINS
  #define NBINS       40     // NUMERO DE BINES
#endif

#ifndef START_BINS
  #define START_BINS  10.0
#endif

#ifndef END_BINS
  #define END_BINS    50000.0
#endif

#ifndef ANGLE_CUT
  #define ANGLE_CUT 60 //M_PI
#endif

#ifdef SKIP
  #ifndef MCUT_MIN
    #define MCUT_MIN 13.5
  #endif
  #ifndef MCUT_MAX
    #define MCUT_MAX 14.0
  #endif
#endif

#ifdef PERCENT
  #ifndef FRACC 
    #define FRACC   0.3           //  EL PORCENTAJE QUE GUARDO
  #endif
#endif

#define NSAMPLES     4  // Number of samples       
#define N_part_types 6  // Number of particle types
#define Ndim         3  // Number of dimension     
#define Nvdir        4  // Number of vdir          

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

extern type_real *tracers;
extern type_real *centros;
extern type_real *vdir_centros;

extern void init_variables(int argc, char **argv);

#endif
