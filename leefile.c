#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <dirent.h> 

#include "variables.h"
#include "cosmoparam.h"
#include "leefile.h"
#include "colores.h"
#include <hdf5.h>
#include <mpi.h>

#define h5_open_file(file_id, name)         H5Fopen(file_id, name, H5P_DEFAULT)
#define h5_open_group(file_id, name)        H5Gopen(file_id, name, H5P_DEFAULT)
#define h5_open_dataset(file_id, name)      H5Dopen(file_id, name, H5P_DEFAULT)
#define h5_open_attribute(parent_id, name)  H5Aopen_name(parent_id, name)

struct io_header header;
struct SnapST snap;

static unsigned long long NumPart_ThisFile(const char *filename, const int type_particle)
{
  unsigned long long npart[N_part_types];

  // Read Header
  hid_t hdf5_file = h5_open_file(filename, H5F_ACC_RDONLY);
  hid_t hdf5_headergrp = h5_open_group(hdf5_file, "/Header");
  hid_t hdf5_attribute = h5_open_attribute(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_ULLONG, &npart);
  H5Aclose(hdf5_attribute);
  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  return npart[type_particle];
}

static void lee(const char *filename, const int MyProc, const int NumProc, const int type_particle, long long *ind)
{
  char key[50];
	long long i, j, pc;
  long long elements;
  herr_t status;
  hid_t hdf5_file;
  hid_t hdf5_typeparticles;
  hid_t hdf5_dataset;
  hid_t hdf5_dataspace;
  hid_t data_type; 
  hsize_t dims[20];
  int rank, type;
  size_t size; 
  void *buff;
  long long numToRead;
  hid_t hdf5_type[6] = {H5T_NATIVE_INT, H5T_NATIVE_LONG, \
                        H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE,\
                        H5T_NATIVE_UINT, H5T_NATIVE_ULONG};
  struct particle_data 
  {
    double      x;
    double      y;
    double      z;
    #ifdef STORE_MASSES
    double      mp;
    #endif
    #ifdef STORE_VELOCITIES
    double      vx;
    double      vy;
    double      vz;
    #endif
    #ifdef STORE_IDS
    unsigned long long id;
    #endif
  } *particles;

  numToRead = (long long)NumPart_ThisFile(filename, type_particle);
  particles = (struct particle_data *) malloc(numToRead*sizeof(struct particle_data));
  if(tracers==NULL)
    tracers = (type_real *) malloc(Ndim*numToRead*sizeof(type_real));
  else
    tracers = (type_real *) realloc(tracers,(long unsigned)Ndim*(*ind + numToRead)*sizeof(type_real));
  assert(tracers!=NULL);

  sprintf(key,"/PartType%d",type_particle);
  hdf5_file = h5_open_file(filename, H5F_ACC_RDONLY);

  for(j=0;j<N_part_types;j++)
  {
    if(j == type_particle)
    { 
      hdf5_typeparticles = h5_open_group(hdf5_file, key);

      hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "Coordinates");
      hdf5_dataspace     = H5Dget_space(hdf5_dataset);
      data_type          = H5Dget_type(hdf5_dataset);
      size               = H5Tget_size(data_type);
      rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
      H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 
      
      elements = 1;
      for(i=0;i<rank;i++)
        elements *= dims[i];

      type = size == 4 ? 2 : 3;
      buff = (void *) malloc(elements*size);
      status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
      
      if(status < 0)
      {
        fprintf(stdout,"Error H5Dread Coordinates\n"); fflush(stdout);
        exit(EXIT_FAILURE);
      }

      for(i=0, pc=0;i<numToRead;i++)
      {
	      switch (type) {
  	      case 2:
            particles[pc].x = ((double)(*(float  *)(buff+size*(3*i  ))))*POSFACTOR;
            particles[pc].y = ((double)(*(float  *)(buff+size*(3*i+1))))*POSFACTOR;
            particles[pc].z = ((double)(*(float  *)(buff+size*(3*i+2))))*POSFACTOR;
		      	break;
  	      case 3: 
            particles[pc].x = ((double)(*(double *)(buff+size*(3*i  ))))*POSFACTOR;
            particles[pc].y = ((double)(*(double *)(buff+size*(3*i+1))))*POSFACTOR;
            particles[pc].z = ((double)(*(double *)(buff+size*(3*i+2))))*POSFACTOR;
            break;
		      default:
		      	fprintf(stdout,"Error! Type Coordinates not correct\n"); fflush(stdout);
  	        exit(EXIT_FAILURE);     
	      }
        pc++;
      }
      
      free(buff);
      H5Tclose(data_type);
      H5Dclose(hdf5_dataset);

#ifdef STORE_VELOCITIES
      hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "Velocities");
      hdf5_dataspace     = H5Dget_space(hdf5_dataset);
      data_type          = H5Dget_type(hdf5_dataset);
      size               = H5Tget_size(data_type);
      rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
      H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 

      elements = 1;
      for(i=0;i<rank;i++)
        elements *= dims[i];

      buff = (void *) malloc(elements*size);
      type = size == 4 ? 2 : 3;
      status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
      if(status < 0)
      {
        fprintf(stdout,"Error H5Dread Velocities\n"); fflush(stdout);
        exit(EXIT_FAILURE);
      }

      for(i=0, pc=0;i<numToRead;i++)
      {
	      switch (type) {
  	      case 2:
            particles[pc].vx = ((double)(*(float  *)(buff+size*(3*i  ))))*VELFACTOR;
            particles[pc].vy = ((double)(*(float  *)(buff+size*(3*i+1))))*VELFACTOR;
            particles[pc].vz = ((double)(*(float  *)(buff+size*(3*i+2))))*VELFACTOR;
		      	break;
  	      case 3:
            particles[pc].vx = ((double)(*(double *)(buff+size*(3*i  ))))*VELFACTOR;
            particles[pc].vy = ((double)(*(double *)(buff+size*(3*i+1))))*VELFACTOR;
            particles[pc].vz = ((double)(*(double *)(buff+size*(3*i+2))))*VELFACTOR;
            break;
		      default:
		      	fprintf(stdout,"Error! Type Velocities not correct\n"); fflush(stdout);
  	        exit(EXIT_FAILURE);     
          	break;
	      }
        pc++;
      }
    
      free(buff);
      H5Tclose(data_type);
      H5Dclose(hdf5_dataset);
#endif
      
#ifdef STORE_MASSES
      hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "Masses");
      hdf5_dataspace     = H5Dget_space(hdf5_dataset);
      data_type          = H5Dget_type(hdf5_dataset);
      size               = H5Tget_size(data_type);
      rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
      H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 

      elements = 1;
      for(i=0;i<rank;i++)
        elements *= dims[i];

      buff = (void *) malloc(elements*size);
      type = size == 4 ? 2 : 3;
      status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
      if(status < 0)
      {
        fprintf(stdout,"Error H5Dread Masses\n"); fflush(stdout);
        exit(EXIT_FAILURE);
      }

      for(i=0, pc=0;i<numToRead;i++)
      {
	      switch (type) {
  	      case 2:
            particles[pc].mp = ((double)(*(float  *)(buff+size*i)));
		      	break;          
  	      case 3:           
            particles[pc].mp = ((double)(*(double *)(buff+size*i)));
            break;
		      default:
		      	fprintf(stdout,"Error! Type Masses not correct\n"); fflush(stdout);
  	        exit(EXIT_FAILURE);     
		      	break;
	      }
        pc++;
      }
    
      free(buff);
      H5Tclose(data_type);
      H5Dclose(hdf5_dataset);
#endif

#ifdef STORE_IDS   
      hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "ParticleIDs");
      hdf5_dataspace     = H5Dget_space(hdf5_dataset);
      data_type          = H5Dget_type(hdf5_dataset);
      size               = H5Tget_size(data_type);
      rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
      H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 

      elements = 1;
      for(i=0;i<rank;i++)
        elements *= dims[i];

      buff = (void *) malloc(elements*size);
      type = size == 4 ? 0 : 1;
      
      // Special case for unsigned integers - treat as unsigned in memory
      H5T_class_t h5class = H5Tget_class(data_type);
      H5T_sign_t  sign  = H5Tget_sign(data_type);
      if(sign==H5T_SGN_NONE && h5class==H5T_INTEGER)
      {
        if(type == 0) 
          type = 4;
        if(type == 1) 
          type = 5;
      }

      // Read dataset
      status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
      if(status < 0)
      {
        fprintf(stdout,"Error H5Dread ParticleIDs\n"); fflush(stdout);
        exit(EXIT_FAILURE);
      }

      for(i=0, pc=0;i<numToRead;i++)
      {
	      switch (type) {
  	      case 0:
            particles[pc].id = ((unsigned long long)(*(int  *)(buff+size*i)));
		      	break;       
  	      case 1:        
            particles[pc].id = ((unsigned long long)(*(long *)(buff+size*i)));
		      	break;       
  	      case 4:        
            particles[pc].id = ((unsigned long long)(*(unsigned int *)(buff+size*i)));
		      	break;       
  	      case 5:        
            particles[pc].id = ((unsigned long long)(*(long unsigned *)(buff+size*i)));
		      	break;
		      default:
		      	fprintf(stdout,"Error! Type ParticleIDs not correct\n"); fflush(stdout);
  	        exit(EXIT_FAILURE);     
		      	break;
	      }
        pc++;
      }
 
      free(buff);
      H5Tclose(data_type);
      H5Dclose(hdf5_dataset);
#endif
      H5Gclose(hdf5_typeparticles);

      if(numToRead != pc)
      {
		     fprintf(stdout,"Error! Bad Read Particles %s - %d - %d\n", filename, pc, numToRead); fflush(stdout);
  	     exit(EXIT_FAILURE);     
      }

    }
  }

  H5Fclose(hdf5_file);

  for(i=0, pc = 0; i<numToRead; i++)
  {
    if(i % NumProc == MyProc)    
    {
#ifdef PERCENT
    if(drand48()<(type_real)FRACC)
    {
#endif
      tracers[((long unsigned)Ndim*(*ind + pc))]     = (type_real)particles[i].x;
      tracers[((long unsigned)Ndim*(*ind + pc)) + 1] = (type_real)particles[i].y;
      tracers[((long unsigned)Ndim*(*ind + pc)) + 2] = (type_real)particles[i].z;
      pc++;
#ifdef PERCENT
    }
#endif
    }
  }

  tracers = (type_real *) realloc(tracers,(long unsigned)Ndim*(*ind + pc)*sizeof(type_real));
  *ind = *ind + pc;
  
  free(particles);

  return;
}

extern void read_gadget(const type_int type_particle [], const int MyProc, const int NumProc)
{
  long long ifile;
  long long ind;
  int index_type;
  char filename[200];
  tracers = NULL;
#ifdef PERCENT
  srand48(80);
#endif

  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++){
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d.hdf5",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s.hdf5",snap.root,snap.name);

    fprintf(stdout,"Leyendo Proc %d - File %s\n",MyProc,filename);
    fflush(stdout);

    for(index_type=0; index_type<NTYPES; index_type++)
      lee(filename,MyProc,NumProc,type_particle[index_type],&ind);
  }

#ifdef PERCENT
  sprintf(message,"Percent %.2f Proc %lu Save %lu representa %.2f %%\n",FRACC,MyProc,(long)ind,100.0f*(type_real)ind/(type_real)cp.npart);
#else
  sprintf(message,"Proc %lu Save %lu representa %.2f %%\n",MyProc,(long)ind,100.0f*(type_real)ind/(type_real)cp.npart);
#endif

  YELLOW(message);
  cp.npart = ind;
  fflush(stdout);

  return;
}

extern void leeheader(const type_int type_particle [])
{
	int i, index_type;
  char filename[200];
  struct io_header header;

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0.hdf5",snap.root,snap.name);
  else
    sprintf(filename,"%s%s.hdf5",snap.root,snap.name);
  
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;

  hdf5_file = h5_open_file(filename, H5F_ACC_RDONLY);
  hdf5_headergrp = h5_open_group(hdf5_file, "/Header");

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &header.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &header.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "Redshift");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "Omega0");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "OmegaLambda");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "HubbleParam");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "Flag_Sfr");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_sfr);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "Flag_Feedback");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_feedback);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = h5_open_attribute(hdf5_headergrp, "Flag_Cooling");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_cooling);
  H5Aclose(hdf5_attribute);

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  cp.npart = 0; 
  cp.Mpart = 0;
  for(index_type = 0; index_type<NTYPES; index_type++)
  {
    i =type_particle[index_type];
    cp.npart += (type_int)header.npartTotal[i] | ((type_int)header.npartTotalHighWord[i] << 32);
    cp.Mpart += header.mass[i];
  }

  // Definicion estructura cosmoparam
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;
  //cp.Mpart    = 3.143E-4*cp.hparam;  /*A-5*/
  //cp.Mpart    = 3.929E-5*cp.hparam;  /*A-4*/
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

  fprintf(stdout,"*********************************** \n");
  fprintf(stdout,"*   Parametros de la simulacion   * \n");
  fprintf(stdout,"*********************************** \n");
  fprintf(stdout,"  Numero de particulas = %llu \n", cp.npart);
  fprintf(stdout,"  Lado del box = %g \n", cp.lbox);
  fprintf(stdout,"  Redshift = %g \n", cp.redshift);
  fprintf(stdout,"  Omega Materia = %g \n", cp.omegam);
  fprintf(stdout,"  Omega Lambda = %g \n", cp.omegal);
  fprintf(stdout,"  Parametro de Hubble = %g \n",cp.hparam);
  fprintf(stdout,"  Masa por particula = %g \n",cp.Mpart);
  fprintf(stdout,"  Softening = %g\n",cp.soft);
  fprintf(stdout,"*********************************** \n");
  fprintf(stdout,"*********************************** \n");
  fflush(stdout);

  return;
}

extern void read_grup_sample(const int MyProc, const int NumProc, const char *filetype)
{
  FILE *pfin;
  char filename[200];
	unsigned long long i, j, k, l, tmp_ngrup;
	unsigned long long npart;
  double mtot;
  double pcm[3];
  double vcm[3];
  double mostbound[3];
  double sig[3];
  double L[3];
  double lambda;
  double m200;
  double r200;
  double v200;
  double mvir;
  double rvir;
  double vvir;
  double vmax;
  double Ep;
  double Ec;
  double aa;
  double bb;
  double cc;
  double evec[3][3];
  double aa_vel;
  double bb_vel;
  double cc_vel;
  double evec_vel[3][3];
  double r;
#ifdef SKIP
	unsigned long long ncontrol = 0;
#endif


  cp.ngrup_sample = 0;
  centros = NULL;
  vdir_centros = NULL;

  k = 0;
  for(l=0; l<NSAMPLES; l++)
  {
#ifdef CM
    sprintf(filename,"../tng_forma/%s_CM_%.2d.bin",filetype,l);
#else
    sprintf(filename,"../tng_forma/%s_%.2d.bin",filetype,l);
#endif
    fprintf(stdout,"Read %d - Filename %s\n",MyProc,filename);

    pfin = fopen(filename,"rb"); 
    fread(&tmp_ngrup,sizeof(unsigned long long),1,pfin);

    fprintf(stdout,"Sample %llu - Ngrups %llu\n",l,tmp_ngrup);
    fflush(stdout);

    cp.ngrup_sample += tmp_ngrup;
    if(centros==NULL)
      centros = (type_real *) malloc(Ndim*cp.ngrup_sample*sizeof(type_real));
    else
      centros = (type_real *) realloc(centros,Ndim*cp.ngrup_sample*sizeof(type_real));
    if(vdir_centros==NULL)
      vdir_centros = (type_real *) malloc(Nvdir*Ndim*cp.ngrup_sample*sizeof(type_real));
    else
      vdir_centros = (type_real *) realloc(vdir_centros,Nvdir*Ndim*cp.ngrup_sample*sizeof(type_real));
    
    assert(centros!=NULL);
    assert(vdir_centros!=NULL);

    for(i=0;i<tmp_ngrup;i++)
    {
	    fread(&npart,sizeof(unsigned long long),1,pfin);
      fread(&mtot,sizeof(double),1,pfin);
      fread(&pcm,sizeof(double),3,pfin);
      fread(&vcm,sizeof(double),3,pfin);
      fread(&mostbound,sizeof(double),3,pfin);
      fread(&sig,sizeof(double),3,pfin);
      fread(&L,sizeof(double),3,pfin);
      fread(&lambda,sizeof(double),1,pfin);
      fread(&m200,sizeof(double),1,pfin);
      fread(&r200,sizeof(double),1,pfin);
      fread(&v200,sizeof(double),1,pfin);
      fread(&mvir,sizeof(double),1,pfin);
      fread(&rvir,sizeof(double),1,pfin);
      fread(&vvir,sizeof(double),1,pfin);
      fread(&vmax,sizeof(double),1,pfin);
      fread(&Ep,sizeof(double),1,pfin);
      fread(&Ec,sizeof(double),1,pfin);
      fread(&aa,sizeof(double),1,pfin);
      fread(&bb,sizeof(double),1,pfin);
      fread(&cc,sizeof(double),1,pfin);
      for(j=0;j<3;j++)
        fread(&evec[j],sizeof(double),3,pfin);
      fread(&aa_vel,sizeof(double),1,pfin);
      fread(&bb_vel,sizeof(double),1,pfin);
      fread(&cc_vel,sizeof(double),1,pfin);
      for(j=0;j<3;j++)
        fread(&evec_vel[j],sizeof(double),3,pfin);

#ifdef SKIP
      mtot = npart > 0 ? 10+log10(mtot) : 0;

      if(mtot < MCUT_MIN || mtot > MCUT_MAX)
      {
        centros[Ndim*k]   = -99;
        centros[Ndim*k+1] = -99;
        centros[Ndim*k+2] = -99;

        for(j=0;j<Nvdir;j++)
        {
          vdir_centros[Ndim*(Nvdir*k+j)]   = -99;
          vdir_centros[Ndim*(Nvdir*k+j)+1] = -99;
          vdir_centros[Ndim*(Nvdir*k+j)+2] = -99;
        }

        k++;
        continue;
      }
#endif

      centros[Ndim*k]   = pcm[0];
      centros[Ndim*k+1] = pcm[1];
      centros[Ndim*k+2] = pcm[2];

      r = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]); 

      vdir_centros[Ndim*Nvdir*k]   = L[0]/r;
      vdir_centros[Ndim*Nvdir*k+1] = L[1]/r;
      vdir_centros[Ndim*Nvdir*k+2] = L[2]/r;

      for(j=1;j<Nvdir;j++)
      {
        vdir_centros[Ndim*(Nvdir*k+j)]   = evec[j-1][0];
        vdir_centros[Ndim*(Nvdir*k+j)+1] = evec[j-1][1];
        vdir_centros[Ndim*(Nvdir*k+j)+2] = evec[j-1][2];
      }

      k++;
	    ncontrol++;
    }

    fclose(pfin);
  }

#ifdef SKIP
  sprintf(message,"To Calculate %llu - %llu // representa %.2f %%\n",ncontrol,k,100.0f*(type_real)ncontrol/(type_real)k);
  YELLOW(message);
#endif

  return;
}
