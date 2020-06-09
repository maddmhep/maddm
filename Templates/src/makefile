include ../Source/make_opts
prog = maddm
fort= $(FC)

libdir = ../lib/
incdir = ../include

flags = $(FFLAGS) -L$(libdir) -I$(incdir)

libs = -lmatrix_elements -ldhelas -lmodel

objs_all=init.o integrate.o interpolate.o maddm.o phasespace.o relic_canon.o relic_coupled.o relic_density.o direct_detection.o directional_detection.o indirect_detect_2to2.o tests.o capture_rate.o alfas_functions.o

#Dependencies
maddm.x: $(incdir)/maddm.inc $(incdir)/maddm_card.inc
direct_detection.o: $(incdir)/maddm.inc $(incdir)/maddm_card.inc
directional_detection.o: $(incdir)/maddm.inc $(incdir)/maddm_card.inc
capture_rate.o: $(incdir)/maddm.inc $(incdir)/capturerate_models.inc
tests.o: $(incdir)/maddm.inc $(incdir)/maddm_card.inc
maddm.o: $(incdir)/maddm.inc $(incdir)/maddm_card.inc

all: maddm

.f.o:
	$(fort) -c $(flags) $*.f

clean:
	-rm -f *.o
	-rm -f ../$(prog).x

maddm: $(objs_all)
	$(fort) $(flags) $(objs_all) $(libs) -o $(prog).x

relic_density: $(objs_all)
	$(fort) $(flags) $(objs_all) $(libs) -o $(prog).x

direct_detection: $(objs_all)
	$(fort) $(flags) $(objs_all) $(libs) -o $(prog).x
