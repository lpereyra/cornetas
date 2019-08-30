### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
#EXTRAS += -DPRECDOUBLE  #Pos and vel in double precision
#EXTRAS += -DLONGIDS     #IDs are long integer
EXTRAS += -DLOG     #BINS LOGARITMICOS
EXTRAS += -DANGLE_CUT=60 #ANGULO DE LA CORNETA
EXTRAS += -DHALO_PARTICULA
EXTRAS += -DPERCENT
EXTRAS += -DFRACC=0.30
EXTRAS += -DSAMPLE_CONTROL  #MUESTRA CONTROL

#CC
CC     := gcc 
CFLAGS := -Wall -O3 -fopenmp -march=native -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm $(GSLL)

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leefile.o variables.o correlacion.o grid.o 

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
