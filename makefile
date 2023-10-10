#
# This is a simple makefile for building Delta_mp2 code
#
#Specify where your MQC_directory is built
MQCDir       = "/home/abovill/MQC_install"

#Compilation with pgfortran
RunF         = pgfortran 
F03Flag      = -mp -i8 -r8 -Mallocatable=03
MQCMODS      = $(MQCDir)/PGI/mod
MQCLIB       = $(MQCDir)/PGI/lib
LIBS         = -llapack -lblas -L$(MQCLIB)

#Compilation with gfortran

#RunF         = gfortran 
#F03Flag      = -fopenmp -fdefault-integer-8 -fdefault-real-8
#MQCMODS      = $(MQCDir)/GNU/mod
#MQCLIB       = $(MQCDir)/GNU/lib
#LIBS         = -llapack -lblas -L$(MQCLIB)
#
#
# The 'all' rule.
#
#
all: mp2.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
%.mod: %.f03
	$(RunF) $(LIBS) $(Prof) -I$(MQCMODS) -c $*.f03

%.o: %.f90
	$(RunF) -I$(MQCMODS) -c $*.f90

%.o: %.f03
	$(RunF) $(F03Flags) -I$(MQCMODS) -c $*.f03

#
# Generic rule for building general executable program (*.exe) from a standard
# f03 source (*.f03) file.
#
%.exe: %.f03 %_mod.f03 
	$(RunF) $(FO3Flag) $(LIBS) -I$(MQCMODS) -o mp2.exe mp2.f03 $(MQCLIB)/libmqc.a
clean:
	rm -rf *.o *.exe *.mod
