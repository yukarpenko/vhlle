
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
EXTRA_FLAGS   = -D SIMPLE  # EoS

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native
LD            = g++
LDFLAGS       = -O3 -march=native

CXXFLAGS     += $(ROOTCFLAGS) $(EXTRA_FLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)

vpath %.cpp src.14.3
objdir     = obj

SRC        = cll.cpp eos.cpp trancoeff.cpp fld.cpp hdo.cpp s95p.cpp ic.cpp main.cpp rmn.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = hlle_visc.14.3
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
