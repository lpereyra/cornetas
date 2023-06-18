#ifndef LEESNAP_H
#define LEESNAP_H

extern void read_grup_sample(const int MyProc, const int NumProc, const char *filetype);
extern void leeheader(const type_int type_particle []);
extern void read_gadget(const type_int type_particle [], const int MyProc, const int NumProc);

/*Input and output files*/
struct SnapST{
  type_int nfiles;
  char root[200], name[200];
  type_int num;
};

struct io_header{
  unsigned int npart[N_part_types];
  double   mass[N_part_types];
  double   time;
  double   redshift;
  unsigned int flag_sfr;
  unsigned int flag_feedback;
  unsigned int npartTotal[N_part_types];
  unsigned int npartTotalHighWord[N_part_types];
  unsigned int flag_cooling;
  unsigned int num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- N_part_types*4- N_part_types*8- 2*8- 2*4- N_part_types*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};

extern struct io_header header;
extern struct SnapST snap;

#endif
