CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
EXTRA_FLAGS   = -D SIMPLE  # EoS

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native
LD            = g++
LDFLAGS       = -O3 -march=native
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS) $(EXTRA_FLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

HYDROO        = cll.o eos.o eo3.o eo1.o eoChiral.o eoHadron.o trancoeff.o fld.o hdo.o s95p.o icurqmd.o ic.o ickw.o main.o rmn.o
 
VPATH=src.14.3
              
TARGET	    = hlle_visc.14.3
#------------------------------------------------------------------------------
$(TARGET):       $(HYDROO)
		$(LD)  $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(HYDROO) $(TARGET)

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $<
	
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<
