### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
#EXTRAS += -DPRECDOUBLE  #Pos and vel in double precision
#EXTRAS += -DLONGIDS     #IDs are long integer
#EXTRAS += -DSTARS
#EXTRAS += -DSTARS_GAS
EXTRAS += -DPERCENT
EXTRAS += -DGROUPS
EXTRAS += -DLOG     #BINS LOGARITMICOS
EXTRAS += -DSKIP    #SALTA LOS HALOS
EXTRAS += -DANGLE_CUT=60 #ANGULO DE LA CORNETA

#CC
CC      := mpicc 
CFLAGS  := -Wall -O3 -march=native -ftree-vectorize -fopenmp -g
HDF5LIB := -lhdf5 #-lz
LIBS    := -lm $(HDF5LIB)

.PHONY : clean 

MAKEFILE := Makefile

OBJS := leefile.o variables.o grid.o correlacion.o 

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := main.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) -c $<

main.x: main.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $^ -o $@ $(LIBS) 

clean:
	rm -rf $(OBJS)
	rm -rf main.o
	rm -rf $(EXEC)
