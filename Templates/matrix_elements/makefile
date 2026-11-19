include ../Source/make_opts
fort= $(FC)

libname = matrix_elements
incdir = ../include

# MadGraph7 matrix elements USE the model_object and aloha_object Fortran modules
moddirs = -I../Source/MODEL -I../Source/DHELAS

flags = $(FFLAGS) -I$(incdir) $(moddirs)

objs = smatrix.o $(patsubst %.f,%.o,$(wildcard *.f))


#MAKEFILE_MADDM --- generates this line when MadDM is run

.f.o:
	$(fort) -c $(flags) $*.f

all: make_lib

clean:
	-rm -f *.o
	-rm -f ../lib/lib$(libname).a

make_lib: $(objs)
	-ar r lib$(libname).a $(objs)
