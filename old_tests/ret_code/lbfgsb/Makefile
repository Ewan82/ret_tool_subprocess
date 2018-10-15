#	Makefile for minimiser LBFGSB
#	txk march 2018

LIBNAME := liblbfgsb-$(FC).a
FC	:= gfortran
FFLAGS	:= -fdefault-real-8

SRCS	:= blas.f  lbfgsb.f  linpack.f  timer.f

OBJS	:= $(SRCS:.f=.o)

$(LIBNAME): $(OBJS)
	ar -r $@ $?

lib: ../liblbfgsb-$(FC).a

clean:
	rm -f *.x $(OBJS) *???~
