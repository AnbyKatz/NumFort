#--------------------------Makefile Template-------------------------
# simply change filename to your filename where ever it appears

# Compiler Flags
F90C = ifort
F90FLAGS = -O2
CFLAGS = -O2 -mkl

# Lapack and numFort libs
LALIBS = -lmkl_lapack95_lp64
NUMLIBS = -I$(NumPath) $(NumPath)kinds.o $(NumPath)lapack.o $(NumPath)numFort.o $(NumPath)minf.o

LIBS = $(NUMLIBS) $(LALIBS)

# Own personal object files
OBJS = 

# Pattern Matching
%.o: %.f90
	$(F90C) $(CFLAGS) -c $< $(LIBS)

#--------------------------Editable MakeFile-------------------------
all: filename

filename: filename.f90 $(OBJS)
	$(F90C) $(CFLAGS) -o $@ $< $(OBJS) $(LIBS) 
clean:
	rm *.o *.mod