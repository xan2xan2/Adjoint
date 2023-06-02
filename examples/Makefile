# Makefile for examples using IBPM library
#
# Clancy Rowley
# Princeton University
#
# $Date$
# $Revision$
# $Author$
# $HeadURL$

include ../config/make.inc

TARGETS = pitching plunging Oseen RigidBodyLoad bin2plt bininfo PO_test_f_NG_8_26 PO_test_f_NG_9_10 PO_test_f_NG_9_10_NO_output PO_test_f_NG_9_11 PO_test_f_NG_9_11_NO_output PO_test_f_NG_9_12 PO_test_f_NG_9_13 PO_test_f_NG_9_11_time PO_test_f_NG_9_11_time_ref PO_test_f_NG_9_17 PO_test_f_NG_9_11_ic_time_ref PO_test_f_NG_9_11_ic_time PO_test_f_NG_22

LDLIBS = -lfftw3 -lm
#MAKEDEPEND = gcc -MM
MAKEDEPEND = mpiCC -MM

LDFLAGS += $(lib_dirs)
#OMPI_LDFLAGS += $(lib_dirs)
CXXFLAGS += -I../src
CXXFLAGS += -I./boost_1_72_0/boost
#OMPI_CXXFLAGS += -I../src

IBPMLIB = libibpm.a

.PHONY: lib clean depend

all: $(TARGETS)

% : %.o lib
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< ../build/$(IBPMLIB) $(LDLIBS)
	#$(OMPI_CXX) $(OPMI_CXXFLAGS) $(OMPI_LDFLAGS) -o $@ $< ../build/$(IBPMLIB) $(LDLIBS)

lib:
	cd ../build && $(MAKE)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $^
	#$(OMPI_CXX) $(OMPI_CXXFLAGS) -c $^

% : %.cc

clean:
	-$(RM) -r *.o $(TARGETS)
	-$(RM) tecplot/*.plt

depend:
	$(MAKEDEPEND) $(CXXFLAGS) *.cc > .depend
	#$(MAKEDEPEND) $(OMPI_CXXFLAGS) *.cc > .depend

-include .depend
