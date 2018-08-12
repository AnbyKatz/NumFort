# Fortran base MakeFile
# Directory to numFort
DIR = /home/anthony/Dropbox/Code/Fortran/f90-toolbox/

F90C = ifort
F90FLAGS = -O2
CFLAGS = -O2 -xHost -mkl

PLPLOT = $$(pkg-config --cflags --libs plplot-fortran)
LIBS = -lmkl_lapack95_lp64 -lplplotfortran -lplplot 
MYLIBS = -I$(DIR)

%.o: %.f90
	$(F90C) $(CFLAGS) -c $<

# Library object files
LIBOBJS = $(DIR)lapack.o\
	$(DIR)numFort.o

# Own personal object files
OBJS = 

all: filename

filename: filename.f90 
	$(F90C) $(CFLAGS) -o $@ $< $(OBJS) $(LIBOBJS) $(LIBS) $(PLPLOT) $(MYLIBS)

clean:
	rm *.o *.mod
