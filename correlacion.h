#ifndef CORRELACION_H
#define CORRELACION_H

extern void correlacion(const int MyProc, const int NumProc, const char *filetype);

struct fcorrST{
  long unsigned ndata_cent;
  long unsigned nrand_cent;
  long unsigned ndata_trac;
  long unsigned nrand_trac;
};

#endif
