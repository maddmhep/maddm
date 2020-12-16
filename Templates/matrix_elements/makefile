include ../Source/make_opts
fort= $(FC)

libname = matrix_elements
incdir = ../include

flags = $(FFLAGS) -I$(incdir)

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
