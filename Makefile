
# The parameter SINGLE_PRECISION makes the code use 'float' instead of
# 'double'.
# The parameter DEBUG enables some additional self-checks of the code.
# The parameter MS_WINDOWS should be used, when the program is compiled for
# MS Windows

.PHONY: all clean cleanall install TAGS
CPPFLAGS = -ggdb -Ioptimisation # -DMS_WINDOWS # -DSINGLE_PRECISION # -DDEBUG
CXXFLAGS = -Wall -O2
LDFLAGS = -L.  # -static -static-libgcc
OBJS1 := coating.o design.o interpolation.o material.o dispersion.o \
        parameters.o reader.o
OBJS2 := analysis.o genetic.o local.o optimisation.o random.o target.o pulse.o
OBJS3 := target.o pulse.o
LDFILES = -lm -lgsl -lgslcblas -lfftw # -lsfftw 

executables = acd acdfield acdinit

all: $(executables)

all_objects = $(OBJS1) $(OBJS2)
$(all_objects): %.o: %.cc %.hh common.hh

acd: $(all_objects) main.cc 
	g++ $(CPPFLAGS) $(CXXFLAGS) $(all_objects) $(LDFLAGS) \
		$(LDFILES) main.cc -o acd

acdfield: $(OBJS1) acdfield.cc
	g++ $(CPPFLAGS) $(CXXFLAGS) $(OBJS1) -lm acdfield.cc -o acdfield

acdinit: $(OBJS1) $(OBJS3) acdinit.cc
	g++ $(CPPFLAGS) $(CXXFLAGS) $(OBJS1) $(OBJS3) acdinit.cc \
	       $(LDFLAGS) $(LDFILES) -o acdinit

install: acd acdfield acdinit
	install acd /usr/local/bin/acd
	install acdfield /usr/local/bin/acdfield
	install acdinit /usr/local/bin/acdinit
	install -d /usr/local/share/acd
	install -m 644 parameters.txt /usr/local/share/acd
	install -m 644 parameters_lite.txt /usr/local/share/acd
	install -d /usr/local/share/acd/materials
	cp -f -R materials/*.dat /usr/local/share/acd/materials/
	chmod 644 /usr/local/share/acd/materials/*.dat

TAGS:
	etags *.hh *.cc

clean:
	rm -f *.o *~ TAGS

cleanall:
	rm -f *.o *~ TAGS $(executables) 
