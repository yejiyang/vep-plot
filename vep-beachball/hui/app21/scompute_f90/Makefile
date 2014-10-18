# ============================================================================
# Name        : Makefile
# Author      : Jiyang Ye
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran

all: src/scompute_f90.f90
	$(FORTRAN_COMPILER) -O2 -g \
		-o bin/scompute_f90 \
		src/scompute_f90.f90

clean:
	rm -f bin/scompute_f90 *.mod
