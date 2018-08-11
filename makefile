# Fortran base MakeFile
# add links to libaries if needed

F90C = ifort
F90FLAGS = -O2
CFLAGS = -O2 -mkl
LIBS = -lmkl_lapack95_lp64

# Path to pre-compiled object files
MYLIBS = -I/home/anthony/Dropbox/Code/Fortran/f90-toolbox/
DIR = /home/anthony/Dropbox/Code/Fortran/f90-toolbox/

%.o: %.f90
	$(F90C) $(CFLAGS) -c $<

# Modules below
LIBOBJS = $(DIR)lapack.o\
	$(DIR)numFort.o

# Add specific object files
OBJS = 

all: fileName

fileName: fileName.f90 $(OBJ)
	$(F90C) $(CFLAGS) -o $@ $< $(OBJS) $(LIBOBJS) $(LIBS) $(MYLIBS)

clean:
	rm *.o *.mod
