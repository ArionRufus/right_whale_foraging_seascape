
# Macros
# the compiler
FC = gfortran -cpp

# GNU ar tool to build static library 
LIBB = ar rc

# flags for debugging or for maximum performance
FCFLAGS = -g -fbounds-check
FCFLAGS = -O2

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/local/netcdf/include

# libraries needed for linking
LDFLAGS += -L/usr/local/netcdf/lib64
LDFLAGS += -lnetcdff -lnetcdf
# variables
OBJECTS = ran_mod.o class_krill.o listKrill_mod.o abstList_mod.o link_mod.o utility_krill_mod.o part_def.o 
UTILITY_OBJ = class_krill.o listKrill_mod.o ran_mod.o part_def.o
LIST_KRILL_OBJ = abstList_mod.o class_krill.o

#MODEL_COUPLED_OBJ = part_def.o class_krill.o listKrill_mod.o utility_krill_mod.o ran_mod.o

vpath %.f90 ../model

# establiishing dependencies and linking files

test: test_2.o  $(OBJECTS)

test_2.o: $(OBJECTS)

utility_krill_mod.o: $(UTILITY_OBJ)

listKrill_mod.o: $(LIST_KRILL_OBJ)

abstList_mod.o: link_mod.o

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<



# utility targets
.PHONY: clean veryclean test1

clean:
	rm -f *.o

veryclean: clean
	rm -f *~ *.a *.mod *.MOD

test1: 
	rm -f *.o *.mod
	make

	./test
