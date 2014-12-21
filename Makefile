PROGRAM = deos
IPOPTLIBDIR = /opt/Ipopt/build/lib
LDFLAGS = -lgfortran -lipopt -lcoinhsl -lcoinlapack -lcoinblas -ldl -lpthread -lgomp
CPPFLAGS = -O3 -pedantic -Wall -std=c++11 -fPIC 
FFLAGS = -O5 -fPIC

CC = g++
FC = gfortran
MAIN = main.cpp
MESHGEN = meshgen

CPPFILES = $(shell find -name '*.cpp' -and -not -name '$(MESHGEN).cpp' -and -not -name '$(MAIN)' -printf '%p ')
FFILES = $(shell find -name '*.f' -printf '%p ')
CPPOBJFILES = $(CPPFILES:%.cpp=%.o)
FOBJFILES = $(FFILES:%.f=%.fo)
CPPDEPFILES = $(CPPFILES:%.cpp=%.d)

all: $(CPPDEPFILES) $(PROGRAM) $(MESHGEN)

$(PROGRAM): $(CPPOBJFILES) $(FOBJFILES) $(MAIN)
	$(CC) -o $(PROGRAM) $(MAIN) $(CPPOBJFILES) $(FOBJFILES) $(LDFLAGS) $(CPPFLAGS)

$(MESHGEN): $(MESHGEN).cpp $(CPPOBJFILES) $(FOBJFILES)
	$(CC) $(CPPFLAGS) -o $(MESHGEN) $(MESHGEN).cpp $(CPPOBJFILES) $(FOBJFILES) $(LDFLAGS) 

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

%.d: %.cpp
	@ $(CC) $(CPPFLAGS) -MM $< -MT $@ | sed -e 's@^\(.*\)\.d:@\1.d \1.o:@' > $@

%.fo: %.f
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: clean 

clean:
	@ rm -f $(CUOBJFILES) $(CPPOBJFILES) $(CPPDEPFILES) $(FOBJFILES) $(PROGRAM) $(MESHGEN)

ifneq "$(MAKECMDGOALS)" "clean"
-include $(CPPDEPFILES) $(CUDEPFILES)
endif

