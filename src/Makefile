ifeq ($(FCOMP),mpifort)
FCFLAGS = -O2 -freal-4-real-8 -march=native -flto
FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
FCBUILD = -Wall -Wextra -pedantic
else ifeq ($(FCOMP),ifort)
FCFLAGS = -O2 -r8 -xHOST -ipo
FCDEBUG = -g -debug -traceback -check all -CB -CU -CA -fpe0 -debug
FCBUILD = -warn all -stand f08
endif

PROGRAM =  ptrace

SRCS =		ran2.f \
			Image_mod.f90 \
			maths.f90 \
			types.f90 \
			main.f90 
                    
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)
debug:	FCFLAGS += $(FCDEBUG)
debug:	$(PROGRAM)
build:	FCFLAGS += $(FCBUILD)
build:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 

%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD ptrace