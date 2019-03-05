## Simple make file, It does NOT check for interdependencies of modules.
## It probably is best to do make clean; make if you change a module.

## Make files are extremely picky about syntax. 
## Indentations have to be done with tabulators!
## Absolutely nothing can follow a continuation backslash.

## List main program module last
## 


SOURCES = precise.f90 \
	constants.f90 \
	input.f90 \
	initial_conditions.f90 \
	print.f90 \
	geo.f90 \
	vorticity.f90 \
	streamfunc.f90 \
	inv.f90 \
	update.f90 \
	SFVE_fdV1-3.f90




## Define name of main program
PROGRAM = test

# Compiler
FF = gfortran

# Delete program
# Linux
RM = rm -f
# DOSe
#RM = del

## Compiler options
# for the Intel Fortran 90 compiler
# CFLAGS = -c -fast -heap-arrays

# for the g95 compiler -g for gdb
CFLAGS = -c -O3 # -02 is safer  # -03 is maybe faster

#Linker Options
# for the Intel Fortran 90 compiler
# LDFLAGS = -fast -heap-arrays

# for the g95 compiler
LDFLAGS = -O2


## Probably no changes necessary below this line
OBJECTS = $(SOURCES:.f90=.o)

all: $(SOURCES) $(PROGRAM)


$(PROGRAM): $(OBJECTS)
	$(FF) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJECTS): %.o: %.f90
	$(FF) $(CFLAGS) $< -o $@ 

clean:
	$(RM) $(OBJECTS) *.mod 

realclean:
	$(RM) $(OBJECTS) *.mod $(PROGRAM)