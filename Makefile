
CXX           = g++
CXXFLAGS      = -Wall -O3 -march=native -flto=auto -fno-math-errno -fopenmp
LD            = g++
LDFLAGS       = -O3 -march=native -flto=auto -fno-math-errno -fopenmp

LIBS          = $(SYSLIBS) -lgsl -lgslcblas

vpath %.cpp src
objdir     = obj

SRC        = cll.cpp eos.cpp eo3.cpp eo1.cpp eoChiral.cpp eoCMF.cpp eoCMFe.cpp eoHadron.cpp eoAZH.cpp eoSmash.cpp \
			 trancoeff.cpp fld.cpp hdo.cpp s95p.cpp icurqmd.cpp ic.cpp ickw.cpp icPartUrqmd.cpp icPartSMASH.cpp \
			 icDynFlu.cpp main.cpp rmn.cpp cornelius.cpp icGlauber.cpp icGubser.cpp icGlissando.cpp icTrento.cpp \
			 icTrento3d.cpp icSuperMC.cpp vtk.cpp icTest.cpp particle.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC))

TARGET	   = hlle_visc
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
