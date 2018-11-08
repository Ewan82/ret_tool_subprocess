#	Makefile for building eispack
#	txk april 2018

FC	:= gfortran
FFLAGS	:= -fdefault-real-8

SRCS	:= eispack.f90 

OBJS	:= $(SRCS:.f90=.o)

LIBNAME = libeis-$(FC).a

.SUFFIXES: .o .f90

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

$(LIBNAME): $(OBJS)
	ar -r $@ $?
#	rm -f $(OBJS)

lib: ../libeis-$(FC).a

clean:
	rm -f *.x *.o *???~
