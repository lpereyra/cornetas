#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leefile.h"
#include "colores.h"

type_real *tracers;
type_real *centros;
type_real *vdir_centros;
struct cosmoparam cp;
char message[200];

extern void init_variables(int argc, char **argv)
{
  FILE *pfin;
  char filename[200];

  RED("Initializing variables...\n");

  sprintf(filename,"%s",argv[1]);
  if(!(pfin=fopen(filename,"r")))
  {
    sprintf(message,"can't open file `%s` \n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%d  \n",&snap.nfiles))
  {
    sprintf(message,"can't read file `%s`\nneed # of snapshots\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%s  \n",snap.root))
  {
    sprintf(message,"can't read file `%s`\nneed snapshots directory\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%s  \n",snap.name))
  {
    sprintf(message,"can't read file `%s`\nneed snapname\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%lf \n",&cp.soft))
  {
    sprintf(message,"can't read file `%s`\nneed softening of simulation\n",filename);RED(message);
    exit(0);
  }

  fclose(pfin);

  /////////////////////////////////////////////////////////////////////////
  char *p = snap.name;
  while (*p) 
  { 
    if(isdigit(*p)) 
    { // Upon finding a digit, ...
        snap.num = strtol(p, &p, 10); // Read a number, ...
    } else { // Otherwise, move on to the next character.
        p++;
    }
  }
  ///////////////////////////////////////////////////////////////////////

  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"Snapnum:                 %d\n",snap.num);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Softening of simulation: %lf \n",cp.soft);BLUE(message);
  BLUE("************* Options ************\n");

  BLUE("********** Makefile Options ***********\n");
  #ifdef PERIODIC
  BLUE("  PERIODIC\n");
  #endif
  #ifdef PRECDOUBLE
  BLUE("  PRECDOUBLE\n");
  #endif
  #ifdef LONGIDS
  BLUE("  LONGIDS\n");
  #endif
  #ifdef POSFACTOR
  sprintf(message,"  POSFACTOR = %f\n",POSFACTOR);BLUE(message);
  #endif
  #ifdef VELFACTOR
  sprintf(message,"  VELFACTOR = %f\n",VELFACTOR);BLUE(message);
  #endif
  #ifdef STORE_VELOCITIES
  BLUE("  STORE_VELOCITIES\n");
  #endif
  #ifdef STORE_IDS
  BLUE("  STORE_IDS\n");
  #endif
  #ifdef STORE_MASSES
  BLUE("  STORE_MASSES\n");
  #endif

  GREEN("END\n");

  return;
}

