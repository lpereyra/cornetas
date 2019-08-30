#ifndef GRID_H
#define GRID_H

#ifndef NGRIDMAX
  #define NGRIDMAX 512
#endif

struct gridst
{
	long nalloc;
	long ngrid;
	type_int nobj;
	type_int *start;
} grid;

extern void grid_init(type_int nobj, type_real *Pos, const type_int cell_size);
extern void grid_free(void);

#endif
