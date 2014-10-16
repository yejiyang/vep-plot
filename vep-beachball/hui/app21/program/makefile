##########################################################################
# Makefile for parallel finite element program generator system          #
##########################################################################
include ddmf.mak
#
FEGEN       =    ..
CC          =    icc
F77         =    ifort
F90         =    ifort
MPIF77      =    mpif90
FLINKER     =    mpif90
OPTFLAGS    =    -O2 -L${FEGEN}/lib
OPTFLAGS    =    -O2 -L${FEGEN}/lib
#OPTFLAGS1   =    -O2 -zero -save -static -L${FEGEN}/lib
#OPTFLAGS1   =    -O2 -zero -save -static -L${FEGEN}/lib
AR          =     ar rv
#
##########################################################################
#
CFLAGS    = -O2 
#CFLAGS    = -O2 -zero -save
FFLAGS    = $(OPTFLAGS)
########## For all kinds of BLAS library  ###########################
###   Fepg blas ###
#LIB_BLAS  = ${FEGEN}/lib/libblas.a
###   Latlas    ###
#LIB_BLAS  = -L/usr/local/atlas/lib -lptf77blas -latlas -lpthread
###   Goto's blas ###
#LIB_BLAS  = -lgotop -lpthread
###  MKL blas   ###
#LIB_BLAS  = -L/opt/intel/mkl/lib/32 -lmkl_p4 -lguide -lpthread
####################################################################
LIB_METIS = ${FEGEN}/lib/libmetis.a 
LIB_FEPGSOLV = ${FEGEN}/lib/libfepgsolv.a
#LIB_LINUX = /usr/lib/gcc-lib/i386-redhat-linux/3.0.4/libg2c.a
#LIB_LINUX = /usr/lib/gcc-lib/i386-redhat-linux/3.2/libg2c.a
#LIB_LINUX = /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/libg2c.a
#LIB_LAPACK =
LIB_MATH  = -lm
LIBS      = $(LIB_BLAS) $(LIB_LAPACK) $(LIB_METIS) $(LIB_FEPGSOLV) $(LIB_LINUX)
# 
COBJ   = io.o partdmesh.o subtime.o
FOBJ   = ${FOBJS} 

default: ddmliq

ddmliq:${COBJ} ${FOBJ}
	$(FLINKER) -o $@ $(FOBJ) $(COBJ) $(LIBS) $(LIB_MATH)
	$(F90) -O2 -o umcazpartition umcazpartition.f io.o partdmesh.o ../lib/libmetis.a
	cp umcazpartition ../compute
	cp ddmliq ../compute/ddmliq21
move: ddmliq 
	mv ddmliq ddmliq.default
.c.o:
	$(CC)  $(CFLAGS) -c $*.c
.F.o:
	$(F90) $(FFLAGS) -c $*.F
.f.o:
	$(MPIF77) $(FFLAGS) -c $*.f
#	$(F90) -O2 -o umcazpartition umcazpartition.o io.o partdmesh.o ../lib/libmetis.a
clall:
	/bin/rm -f   *.o ddm.o ddmliq ddmliq.default
clean:
	/bin/rm -f   *.o
#
#========================================================================
#
