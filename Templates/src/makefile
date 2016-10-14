prog = maddm
fort= gfortran

libdir = ../lib/
incdir = ../include

flags = -O -fno-align-commons -ffixed-line-length-132 -L$(libdir) -I$(incdir)

libs = -lmatrix_elements -ldhelas -lmodel

objs_all=init.o integrate.o interpolate.o maddm.o phasespace.o relic_canon.o relic_coupled.o relic_density.o direct_detection.o directional_detection.o tests.o

#Dependencies
init.o: $(incdir)/maddm.inc
integrate.o: $(incdir)/maddm.inc
interpolate.o: $(incdir)/maddm.inc
maddm.x: $(incdir)/maddm.inc $(incdir)/maddm_card.inc
phasespace.o: $(incdir)/maddm.inc
relic_canon.o: $(incdir)/maddm.inc
relic_coupled.o: $(incdir)/maddm.inc
relic_density.o: $(incdir)/maddm.inc
direct_detection.o: $(incdir)/maddm.inc
directional_detection.o: $(incdir)/maddm.inc
tests.o: $(incdir)/maddm.inc

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
