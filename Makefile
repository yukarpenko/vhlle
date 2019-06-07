# hs == hadron sampler (particlization)
# jt == parton cascade
OBJ = obj/
hDir = src/
jtDir = JT/
hsDir = HS/

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

CXX           = g++
CXXFLAGS      = -std=c++11 -Wall -fPIC -O3 -march=native -Isrc/ -IJT/ -I../pythia8235/include/
LD            = g++
LDFLAGS       = -O3 -march=native

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) ../pythia8235/lib/libpythia8.a -ldl

#vpath %.cpp src JT
#objdir     = obj

hydro_src = $(hDir)cll.cpp $(hDir)eos.cpp $(hDir)eo3.cpp $(hDir)eo1.cpp \
       $(hDir)eoChiral.cpp $(hDir)eoHadron.cpp $(hDir)eoAZH.cpp $(hDir)trancoeff.cpp \
       $(hDir)fld.cpp $(hDir)hdo.cpp $(hDir)s95p.cpp $(hDir)icurqmd.cpp $(hDir)ic.cpp \
       $(hDir)ickw.cpp $(hDir)icPartUrqmd.cpp $(hDir)main.cpp $(hDir)rmn.cpp \
       $(hDir)cornelius.cpp $(hDir)icGlauber.cpp $(hDir)icGubser.cpp \
       $(hDir)icGlissando.cpp $(hDir)icEpos.cpp $(hDir)smearingKernel.cpp
jt_src = $(jtDir)params.cpp $(jtDir)integral.cpp $(jtDir)W.cpp $(jtDir)S.cpp \
         $(jtDir)cascade.cpp $(jtDir)interfacePythia.cpp
hs_src = $(hsDir)cascade.cpp  $(hsDir)DecayChannel.cpp  $(hsDir)particle.cpp \
         $(hsDir)tree.cpp  $(hsDir)DatabasePDG2.cpp  $(hsDir)gen.cpp \
         $(hsDir)params.cpp  $(hsDir)ParticlePDG2.cpp  $(hsDir)UKUtility.cpp

V : dirs  hlle_visc
#------------------------------------------------------------------------------
hlle_visc:  $(hydro_src:%.cpp=$(OBJ)%.o) $(jt_src:%.cpp=$(OBJ)%.o) $(hs_src:%.cpp=$(OBJ)%.o)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
$(hydro_src:%.cpp=$(OBJ)%.o) : $(OBJ)%.o : ./%.cpp
	$(CXX) $(CXXFLAGS)  -c $< -o $@
$(jt_src:%.cpp=$(OBJ)%.o) : $(OBJ)%.o : ./%.cpp
	$(CXX) $(CXXFLAGS)  -c $< -o $@
$(hs_src:%.cpp=$(OBJ)%.o) : $(OBJ)%.o : ./%.cpp
	$(CXX) $(CXXFLAGS)  -c $< -o $@
dirs:
	@if [ ! -d $(OBJ)  ]  ;then                       \
	set -x; mkdir $(OBJ); set +x;                  \
	fi
	@if [ ! -d $(OBJ)$(hDir)  ]  ;then                       \
	set -x; mkdir  $(OBJ)$(hDir); set +x;                  \
	fi
	@if [ ! -d $(OBJ)$(jtDir)  ]  ;then                       \
	set -x; mkdir  $(OBJ)$(jtDir); set +x;                  \
	fi
	@if [ ! -d $(OBJ)$(hsDir)  ]  ;then                       \
	set -x; mkdir  $(OBJ)$(hsDir); set +x;                  \
	fi
clean:
		@rm -f $(OBJ)$(hDir)/*.o $(OBJ)$(jtDir)/*.o $(OBJ)$(hsDir)/*.o hlle_visc
