#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h> 
#include <dirent.h> 
#include <math.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leefile.h"
#include "colores.h"

extern void leeheader(void)
{

  FILE *pf;
  type_int d1,d2;
  char filename[200];

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fclose(pf);

  /* Definicion estructura cosmoparam */
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize*POSFACTOR;
  cp.npart     = header.npartTotal[1];
  cp.Mpart     = header.mass[1];
  //cp.Mpart    = 3.143E-4*cp.hparam;  /*A-5*/
  //cp.Mpart    = 3.929E-5*cp.hparam;  /*A-4*/
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

  printf("*********************************** \n");
  printf("*   Parametros de la simulacion   * \n");
  printf("*********************************** \n");
  printf("  Numero de particulas = %u \n", cp.npart);
  printf("  Lado del box = %g \n", cp.lbox);
  printf("  Redshift = %g \n", cp.redshift);
  printf("  Omega Materia = %g \n", cp.omegam);
  printf("  Omega Lambda = %g \n", cp.omegal);
  printf("  Parametro de Hubble = %g \n",cp.hparam);
  printf("  Masa por particula = %g \n",cp.Mpart);
  printf("  Softening = %g\n",cp.soft);
  printf("*********************************** \n");
  printf("*********************************** \n");
}


static int compare_str(const void* a, const void* b) 
{ 
  return strcmp(*(const char**)a, *(const char**)b); 
}

static char *set_filename(const char *folder, const char *search, const type_int FILE_SAMPLE) 
{ 
    int i;
    struct dirent *de;  // Pointer for directory entry 
	  const char delim[] = "_";
    char *str, *tmp, *ptr;
  	char **char_arr;
  	int  *size_arr;
    DIR *dr;
 
    // opendir() returns a pointer of DIR type.  
    size_arr = (int *) calloc(NSAMPLES,sizeof(int));
  
    dr = opendir(folder);
    assert(dr != NULL);

    i = 0;
    while((de = readdir(dr)) != NULL)
    {
      if(strstr(de->d_name, search) != NULL)
      {
        i++;
      }
    }

    assert(i==NSAMPLES);

    dr = opendir(folder);
    assert(dr != NULL);

    i = 0;
    while((de = readdir(dr)) != NULL)
    {
      if(strstr(de->d_name, search) != NULL)
      {
        size_arr[i] = strlen(de->d_name);
        i++;
      }
    }

    char_arr = (char **) malloc(NSAMPLES*sizeof(char *));
    for(i=0;i<NSAMPLES;i++)
      char_arr[i] = (char *) malloc(size_arr[i]*sizeof(char));

    dr = opendir(folder);
    assert(dr != NULL);

    i = 0;
    while((de = readdir(dr)) != NULL)
    {
      if(strstr(de->d_name, search) != NULL)
      {
        strcpy(char_arr[i],de->d_name);
        i++;
      }
    }

    qsort(char_arr, NSAMPLES, sizeof(const char *), compare_str);

    str = malloc(char_arr[FILE_SAMPLE]);
    strcpy(str, char_arr[FILE_SAMPLE]);
  	ptr = strtok(str, delim);

    i = 0;
  	while(ptr != NULL)
	  {
      if(i==1)
      {
        tmp = malloc(strlen(ptr));
        strcpy(tmp, ptr);
        i++;

      }else if(i>1){

        snap.sample_info = malloc(strlen(tmp) + strlen(delim) + strlen(ptr) + 1);
        strcpy(snap.sample_info, tmp);
        strcat(snap.sample_info, delim);
        strcat(snap.sample_info, ptr);

        tmp = malloc(strlen(snap.sample_info));
        strcpy(tmp, snap.sample_info);
      }

      if(strstr(ptr, "sample") != NULL)
      {
        i++;
      }

  		ptr = strtok(NULL, delim);
	  }

    tmp = snap.sample_info;
    tmp[strlen(tmp) - 4] = 0; // le quito la extension ".bin"
    strcpy(snap.sample_info, tmp);

    tmp = malloc(strlen(folder) + strlen(char_arr[FILE_SAMPLE]));
    strcpy(tmp, folder);
    strcat(tmp, char_arr[FILE_SAMPLE]);

    return tmp;
} 

extern void read_grup_sample(type_real *fof, const type_int FILE_SAMPLE)
{
  char  *filename;
  type_int  i, j;
  FILE  *pfin;
  type_int save;
  type_int id;
  type_int NumPart;
  type_real Mass;
  type_real aa;
  type_real bb;
  type_real cc;
  type_real ParP;
  type_real vcm[3];
  type_real L[3];
  #ifdef SAMPLE_CONTROL
    type_real vdir[3];
    filename = set_filename("../samples/","halos_control", FILE_SAMPLE); 
  #else
    filename = set_filename("../samples/","halos_in_fil", FILE_SAMPLE); 
  #endif
  type_real r;

  pfin = fopen(filename,"rb"); 
  fprintf(stdout,"Read %s\n",filename);

  fread(&cp.ngrup_sample,sizeof(type_int),1,pfin);        
  centros = (type_real *)  malloc(Ndim*cp.ngrup_sample*sizeof(type_real));
  vdir_centros = (type_real *) malloc(Nvdir*Ndim*cp.ngrup_sample*sizeof(type_real));

  fprintf(stdout,"Grupos in sample %d\n",cp.ngrup_sample);
  fflush(stdout);

  for(i=0;i<cp.ngrup_sample;i++)
  {
    fread(&save,sizeof(type_int),1,pfin);
    fread(&id,sizeof(type_int),1,pfin);
    fread(&centros[Ndim*i],sizeof(type_real),1,pfin);
    fread(&centros[Ndim*i+1],sizeof(type_real),1,pfin);
    fread(&centros[Ndim*i+2],sizeof(type_real),1,pfin);
    fread(&NumPart,sizeof(type_int),1,pfin);
    fread(&Mass,sizeof(type_real),1,pfin);
    fread(&aa,sizeof(type_real),1,pfin);
    fread(&bb,sizeof(type_real),1,pfin);
    fread(&cc,sizeof(type_real),1,pfin);
    fread(&ParP,sizeof(type_real),1,pfin);
    fread(&vcm[0],sizeof(type_real),3,pfin);
    fread(&L[0],sizeof(type_real),3,pfin);

    r = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]); 

    vdir_centros[Ndim*Nvdir*i]   = L[0]/r;
    vdir_centros[Ndim*Nvdir*i+1] = L[1]/r;
    vdir_centros[Ndim*Nvdir*i+2] = L[2]/r;

    for(j=1;j<Nvdir;j++)
    {
      fread(&vdir_centros[Ndim*(Nvdir*i+j)]  ,sizeof(type_real),1,pfin);
      fread(&vdir_centros[Ndim*(Nvdir*i+j)+1],sizeof(type_real),1,pfin);
      fread(&vdir_centros[Ndim*(Nvdir*i+j)+2],sizeof(type_real),1,pfin);
    }

    #ifdef SAMPLE_CONTROL
    fread(&vdir[0],sizeof(type_real),1,pfin);
    fread(&vdir[1],sizeof(type_real),1,pfin);
    fread(&vdir[2],sizeof(type_real),1,pfin);
    #endif
  }

  fclose(pfin);

  return;

}

#ifdef HALO_PARTICULA

static void lee(const char *filename, type_int *ind){
  FILE *pf;
  type_int d1, d2;
  type_int k, pc, n;
  type_real r[3];
  #ifdef STORE_VELOCITIES
  type_real v[3];
  #endif
  #ifdef STORE_IDS
  type_int id;
  #endif
#ifdef PERCENT
  type_real *Q;
  //struct particle_data *Q;
#endif

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  
#ifdef PERCENT
  //Q = (struct particle_data *) malloc(header.npart[1]*sizeof(struct particle_data));
  Q = (type_real *) malloc(Ndim*header.npart[1]*sizeof(type_real));
  assert(Q != NULL);
#endif

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&r[0], sizeof(type_real), 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
#ifdef PERCENT
        Q[Ndim*pc]   = r[0]*POSFACTOR;
        Q[Ndim*pc+1] = r[1]*POSFACTOR;
        Q[Ndim*pc+2] = r[2]*POSFACTOR;
        //Q[pc].Pos[0] = r[0]*POSFACTOR;
        //Q[pc].Pos[1] = r[1]*POSFACTOR;
        //Q[pc].Pos[2] = r[2]*POSFACTOR;
#else
        P[*ind+pc].Pos[0] = r[0]*POSFACTOR;
        P[*ind+pc].Pos[1] = r[1]*POSFACTOR;
        P[*ind+pc].Pos[2] = r[2]*POSFACTOR;
#endif
        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_VELOCITIES
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&v[0], sizeof(type_real), 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
#ifdef PERCENT
        Q[pc].Vel[0] = v[0]*VELFACTOR;
        Q[pc].Vel[1] = v[1]*VELFACTOR;
        Q[pc].Vel[2] = v[2]*VELFACTOR;
#else
        P[*ind+pc].Vel[0] = v[0]*VELFACTOR;
        P[*ind+pc].Vel[1] = v[1]*VELFACTOR;
        P[*ind+pc].Vel[2] = v[2]*VELFACTOR;
#endif
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_IDS
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&id, sizeof(type_int), 1, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
#ifdef PERCENT
        Q[pc].id = id;
#else
        P[*ind+pc].id = id;
#endif
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

#ifdef PERCENT

  //if(P==NULL)
  //  P = (struct particle_data *) malloc(header.npart[1]*sizeof(struct particle_data));
  //else
  //  P = (struct particle_data *) realloc(P,(*ind + header.npart[1])*sizeof(struct particle_data));
  //assert(P!=NULL);

  if(tracers==NULL)
    tracers = (type_real *) malloc(Ndim*header.npart[1]*sizeof(type_real));
  else
    tracers = (type_real *) realloc(tracers,Ndim*(*ind + header.npart[1])*sizeof(type_real));
  assert(tracers!=NULL);

  pc = 0;
  for(n=0; n<header.npart[1]; n++)
  {
    if(drand48()<FRACC)
    {
      //P[*ind + pc] = Q[n];
      tracers[Ndim*(*ind + pc)]   = Q[Ndim*n];
      tracers[Ndim*(*ind + pc)+1] = Q[Ndim*n+1];
      tracers[Ndim*(*ind + pc)+2] = Q[Ndim*n+2];
      pc++;
    }
  }

  free(Q);

  //P = (struct particle_data *) realloc(P,(*ind + pc)*sizeof(struct particle_data));
  tracers = (type_real *) realloc(tracers,Ndim*(*ind + pc)*sizeof(type_real));

#endif

  *ind += pc;
  
  fclose(pf);
}

extern void read_gadget(void){
  char filename[200];
  type_int ifile, ind;
#ifdef PERCENT

    srand48(80);

    tracers = NULL;
    //P = NULL;

#else
    size_t total_memory;

    /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
    total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
    printf("Allocating %.5zu Gb for %u particles\n",total_memory,cp.npart);
    P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
    assert(P != NULL);
 
#endif

  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++){
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,&ind);
  }
  
#ifdef PERCENT

  sprintf(message,"De %lu me guardo %lu, representa %.2f %%\n",(long)cp.npart,(long)ind,100.0f*(type_real)ind/(type_real)cp.npart);
  YELLOW(message);
  fflush(stdout);
  cp.npart = ind;

#endif

  fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);
}

#else

extern void read_grup_fof(type_real *fof)
{
  char  filename[200];
  type_int  i, k;
  FILE  *pfin;
  type_int save;
  type_int id;
  type_int NumPart;
 
  #ifdef MCRITIC
    sprintf(filename,"../%.2d_%.2f_centros_cut_%.2f.bin",snap.num,fof[1],m_critica);
  #else
    sprintf(filename,"../../%.2d_%.2f_centros.bin",snap.num,fof[1]);
  #endif

  pfin = fopen(filename,"rb"); 

  fread(&cp.ngrup,sizeof(type_int),1,pfin);

  fprintf(stdout,"Grupos %d\n",cp.ngrup);
  fflush(stdout);

  srand48(80);

  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

  for(i=0;i<cp.ngrup;i++)
  {
    fread(&save,sizeof(type_int),1,pfin);
    fread(&id,sizeof(type_int),1,pfin);
    fread(&tracers[Ndim*i],  sizeof(type_real),1,pfin);
    fread(&tracers[Ndim*i+1],sizeof(type_real),1,pfin);
    fread(&tracers[Ndim*i+2],sizeof(type_real),1,pfin);
    fread(&NumPart,sizeof(type_int),1,pfin);
  }

  fclose(pfin);

  return;

}

#endif
