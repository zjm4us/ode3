ODELIB = .

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)

CXXFLAGS  += $(ROOTCFLAGS) -I$(ODELIB) -Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(ODELIB) -L$(ODELIB)
GXX        = g++

EXES = baseball1 baseball2 RKnDemo RKnStep RKnTest projGSL

all: $(EXES)

RKnDemo: RKnDemo.cpp RKn.cpp $(ODELIB)/RKn.hpp
	$(GXX) $(CXXFLAGS) -o $@ RKnDemo.cpp RKn.cpp $(LDFLAGS)

RKnStep: RKnStep.cpp RKn.cpp $(ODELIB)/RKn.hpp
	$(GXX) $(CXXFLAGS) -o $@ RKnStep.cpp RKn.cpp $(LDFLAGS)

RKnTest: RKnTest.cpp RKn.cpp $(ODELIB)/RKn.hpp
	$(GXX) $(CXXFLAGS) -o $@ RKnTest.cpp RKn.cpp $(LDFLAGS)

baseball1: baseball1.cpp RKn.cpp $(ODELIB)/RKn.hpp
	$(GXX) $(CXXFLAGS) -o $@ baseball1.cpp RKn.cpp $(LDFLAGS)

baseball2: baseball2.cpp RKn.cpp $(ODELIB)/RKn.hpp
	$(GXX) $(CXXFLAGS) -o $@ baseball2.cpp RKn.cpp $(LDFLAGS)

projGSL: projGSL.cpp
	$(GXX) $(CXXFLAGS) -o $@ projGSL.cpp $(LDFLAGS) $(GSLFLAGS)

vterm: vterm.cpp RKn.cpp $(ODELIB)/RKn.hpp
	$(GXX) $(CXXFLAGS) -o $@ vterm.cpp RKn.cpp $(LDFLAGS)

vterm_energy: vterm_energy.cpp RKn.cpp RKn.hpp
	$(CXX) $(CXXFLAGS) -o $@ vterm_energy.cpp RKn.cpp $(LDFLAGS)

vterm_air: vterm_air.cpp RKn.cpp RKn.hpp
	$(CXX) $(CXXFLAGS) -o $@ vterm_air.cpp RKn.cpp $(LDFLAGS)

vterm_vs_mass: vterm_vs_mass.cpp
	$(CXX) $(CXXFLAGS) `root-config --cflags --libs` -o vterm_vs_mass vterm_vs_mass.cpp



clean:
	rm -f *.o *.so *.d *~ $(EXES) RKnDemo.root Projectile.png

