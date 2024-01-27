FC = gfortran
RM = rm -f

FFLAGS = -O3

EXE=modir_dfmo
OBJS=halton.o sobol.o qsortd.o latin_random.o
OBJS_DFMO=main_DFMO.o subroutines_DFMO.o DFMO_strong_sottomesso.o
OBJS_MODIR=main_MODIRECT.o  mainSUB.o  struttura_dati.o  subroutines_DIRECT.o
MODS=modules_DFMO.o modules_DIRECT.o 

all:  $(MODS) main.o $(OBJS) $(OBJS_DFMO) subroutines_RACCORDO.o $(OBJS_MODIR) problem.o
	$(FC) -o $(EXE) $(MODS) main.o $(OBJS) $(OBJS_DFMO) subroutines_RACCORDO.o $(OBJS_MODIR) problem.o 
	@echo "\nDone compiling $(EXE)\n"

.SUFFIXES: .f90 .o
	.f   .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90

.f.o:   $* ; $(FC) $(FFLAGS) -c $*.f

clean: 
	$(RM) *.o
	$(RM) *.mod
	$(RM) $(EXE)
	$(RM) best.txt
	$(RM) out_nondom.txt
	$(RM) fort.*
