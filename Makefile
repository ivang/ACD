
# The parameter SINGLE_PRECISION makes the code use 'float' instead of
# 'double'.
# The parameter DEBUG enables some additional self-checks of the code.
# The parameter MS_WINDOWS should be used, when the program is compiled for
# MS Windows

.PHONY: all clean cleanall install uninstall TAGS

CC = g++
CPPFLAGS = -ggdb -Ioptimisation # -DMS_WINDOWS # -DSINGLE_PRECISION # -DDEBUG
CXXFLAGS = -Wall #-O2
LDFLAGS = -L.  # -static -static-libgcc
OBJS1 := coating.o design.o interpolation.o material.o dispersion.o \
        parameters.o reader.o
OBJS2 := analysis.o genetic.o local.o optimisation.o random.o target.o pulse.o
OBJS3 := target.o pulse.o
LDFILES = -lm -lgsl -lgslcblas -lfftw # -lsfftw 

INSTALL = install
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA = $(INSTALL) -m 644

RM = rm -f
TAGGER = ctags

executables = acd acdfield acdinit

prefix = /usr/local
bindir = $(prefix)/bin
datadir = $(prefix)/share/acd

all: $(executables)

all_objects = $(OBJS1) $(OBJS2)
$(all_objects): %.o: %.cc %.hh common.hh

acd: $(all_objects) main.cc 
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(all_objects) $(LDFLAGS) \
		$(LDFILES) main.cc -o acd

acdfield: $(OBJS1) acdfield.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(OBJS1) $(LDFILES) acdfield.cc -o acdfield

acdinit: $(OBJS1) $(OBJS3) acdinit.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(OBJS1) $(OBJS3) acdinit.cc \
	       $(LDFLAGS) $(LDFILES) -o acdinit

install: acd acdfield acdinit
	$(INSTALL) acd acdfield acdinit $(bindir)
	$(INSTALL) -d $(datadir)
	$(INSTALL_DATA) parameters.txt $(datadir)
	$(INSTALL_DATA) parameters_lite.txt $(datadir)
	$(INSTALL) -d $(datadir)/materials
	$(INSTALL_DATA) materials/*.dat $(datadir)/materials

uninstall:
	$(RM) $(bindir)/acd $(bindir)/acdfield $(bindir)/acdinit
	$(RM) $(datadir)/parameters.txt $(datadir)/parameters_lite.txt
	$(RM) $(datadir)/materials/*.dat
	rmdir $(datadir)/materials/
	rmdir $(datadir)

TAGS:
	$(TAGGER) *.hh *.cc

clean:
	rm -f *.o *~ TAGS tags

cleanall:
	rm -f *.o *~ TAGS tags  $(executables)


