prog = maddm
fort= gfortran

libdir = ../lib/
incdir = ../include

flags = -O -fno-align-commons -ffixed-line-length-132 -L$(libdir) -I$(incdir)

libs = -lmatrix_elements -ldhelas -lmodel

objs_all=init.o integrate.o interpolate.o maddm.o phasespace.o relic_canon.o relic_coupled.o relic_density.o direct_detection.o directional_detection.o tests.o

all: clean maddm

.f.o:
	$(fort) -c $(flags) $*.f

clean:
	-rm -f *.o
	-rm -f ../$(prog).x

maddm: $(objs_all)
	$(fort) $(flags) $(objs_all) $(libs) -o $(prog).x

relic_density: clean $(objs_all)
	$(fort) $(flags) $(objs_all) $(libs) -o $(prog).x

direct_detection:  clean $(objs_all)
	$(fort) $(flags) $(objs_all) $(libs) -o $(prog).x
