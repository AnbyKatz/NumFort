# Directory to numFort
# DO NOT USE ~ FOR THE HOME DIRECTORY IN THIS PATH
DIR = /path/to/numFort/

F90C = ifort
F90FLAGS = -O2
CFLAGS = -O2 -mkl

LIBS = -lmkl_lapack95_lp64
MYLIBS = -I$(DIR)
# Comment out PLPLOT if it wasn't installed and all of its appearance's
PLPLOT = -lplplotfortran -lplplot $$(pkg-config --cflags --libs plplot-fortran)

%.o: %.f90
	$(F90C) $(CFLAGS) -c $<

# Library object files
LIBOBJS = $(DIR)kinds.o\
	$(DIR)lapack.o\
	$(DIR)numFort.o\
	$(DIR)PLplots.o

# Own personal object files
OBJS = 

all: filename

filename: filename.f90 $(OBJS)
	$(F90C) $(CFLAGS) -o $@ $< $(OBJS) $(LIBOBJS) $(LIBS) $(MYLIBS) $(PLPLOT)

# Dependencies

clean:
	rm *.o *.mod
