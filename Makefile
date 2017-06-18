
ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs`

CXXFLAGS += -I. -Wall


%.o : %.c
	$(RM) $@
	$(CC) -c $(CFLAGS) -o $@ $<

%.o : %.cc
	$(RM) $@
	$(CXX) -c $(CXXFLAGS) -o $@ $*.cc


OBJS    = EarthDensity.o BargerPropagator.o mosc.o mosc3.o


LIBBASE   = ThreeProb
VER       = 2.10
TAG       =
LIBALIAS  = $(LIBBASE)$(TAG)
LIBNAME   = $(LIBALIAS)_$(VER)

lib3p     = lib$(LIBNAME).a


targets = $(lib3p) probRoot probLinear probAnalytic


$(lib3p) : $(OBJS)
	$(RM) $@
	ar clq $@ $(OBJS)
	ranlib $@


probRoot: probRoot.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probRoot.o
probRoot.o:
	$(CXX) -o probRoot.o $(ROOTCFLAGS) $(CXXFLAGS) -c probRoot.cc



probLinear: probLinear.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probLinear.o
probLinear.o:
	$(CXX) -o probLinear.o $(ROOTCFLAGS) $(CXXFLAGS) -c probLinear.cc



oscillateOne: oscillateOne.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: oscillateOne.o
oscillateOne.o:
	$(CXX) -o oscillateOne.o $(ROOTCFLAGS) $(CXXFLAGS) -c oscillateOne.cc

os1: os1.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: os1.o
os1.o:
	$(CXX) -o os1.o $(ROOTCFLAGS) $(CXXFLAGS) -c os1.cc

# MY (ANDREW YONG) OWN HEADER
os1_yong: os1_yong.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: os1_yong.o
os1_yong.o:
		$(CXX) -o os1_yong.o $(ROOTCFLAGS) $(CXXFLAGS) -c os1_yong.cc


os1h_yong: os1h_yong.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: os1h_yong.o
os1h_yong.o:
		$(CXX) -o os1h_yong.o $(ROOTCFLAGS) $(CXXFLAGS) -c os1h_yong.cc


fubar: fubar.o Nu_Fitter.o Appearance.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: fubar.o
fubar.o:
		$(CXX) -o fubar.o $(ROOTCFLAGS) $(CXXFLAGS) -c fubar.cpp 

.PHONY: Nu_Fitter.o
Nu_Fitter.o:
		$(CXX) -o Nu_Fitter.o $(ROOTCFLAGS) $(CXXFLAGS) -c Nu_Fitter.cpp

.PHONY: Appearance.o
Appearance.o:
		$(CXX) -o fubar.o $(ROOTCFLAGS) $(CXXFLAGS) -c Appearance.cpp



# END OF ANDREW YONG HEADER

CreateHistsToFile: CreateHistsToFile.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: CreateHistsToFile.o
CreateHistsToFile.o:
	$(CXX) -o CreateHistsToFile.o $(ROOTCFLAGS) $(CXXFLAGS) -c CreateHistsToFile.cc




probAnalytic: probAnalytic.o $(lib3p)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probAnalytic.o
probAnalytic.o:
	$(CXX) -o probAnalytic.o $(ROOTCFLAGS) $(CXXFLAGS) -c probAnalytic.cc



.PHONY: all
all: $(targets)


.PHONY: clean
clean:
	$(RM) $(targets) *.o

emptyrule:: $(lib3p)
