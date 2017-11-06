export

SHELL = /bin/sh

ARCH = x86_64
OS   = linux-gnu

HOST = castor0
SITE = cscs.ch

srcdir = .


prefix = .
exec_prefix = build/x86_64-unknown-linux-gnu

bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib

NETCDFROOT     = /apps/eiger/netcdf/4.13/gnu_434
NETCDFLIBPATH  = $(NETCDFROOT)/lib
NETCDF_LIB     = -L$(NETCDFLIBPATH) -lnetcdff -lnetcdf
NETCDF_INCLUDE = -I$(NETCDFROOT)/include

HDF5ROOT       = 
HDF5_LIB       =
HDF5_INCLUDE   =

SZIPROOT       = 
SZIP_LIB       = 
SZIP_INCLUDE   = 

ZLIBROOT       = 
ZLIB_LIB       = 
ZLIB_INCLUDE   = 

MPIROOT        = 
MPI_LIB        = 
MPI_INCLUDE    = 

LAPACKROOT     = 
LAPACK_LIB     = 


PROFILE_LIB     = 
PROFILE_INCLUDE = 

OPENCL_LIB     = -lOpenCL -lstdc++
OPENCL_INCLUDE = -I/apps/eiger/Cuda-4.0/cuda/include

LIBS           = -L../lib -lsupport $(LAPACK_LIB) $(NETCDF_LIB) $(HDF5_LIB) $(SZIP_LIB) $(ZLIB_LIB) $(MPI_LIB) $(PROFILE_LIB) $(OPENCL_LIB)

INCLUDE        = -I../include $(NETCDF_INCLUDE) $(HDF5_INCLUDE) $(SZIP_INCLUDE) $(ZLIB_INCLUDE) $(MPI_INCLUDE) $(PROFILE_INCLUDE) $(OPENCL_INCLUDE)
INCLUDES       = $(INCLUDE) 

AS             = as

CC             = gcc
CFLAGS         = $(INCLUDE) -march=native -O -DpgiFortran -DHAVE_LIBNETCDF -DHAVE_CF_INTERFACE -D__ICON_TESTBED_BASE__ -D__ICON__ -D__ICON_TESTBED__ -DNOMPI -DTURBO_MODE
FC             = gfortran
FFLAGS         = $(INCLUDES) -J../module -I../module -march=native -O -ffast-math -xf95-cpp-input -std=f2003 -fmodule-private -fimplicit-none -fmax-identifier-length=31 -ffree-line-length-99 -Wall -finit-real=nan -D__LOOP_EXCHANGE -D__ICON_TESTBED_BASE__ -D__ICON__ -D__ICON_TESTBED__ -DNOMPI
F77            = gfortran
F77FLAGS       = -march=native -ffast-math

AR             = ar
ARFLAGS        = crv

LDFLAGS        = -march=native -ffast-math -O

SRCDIRS        =  support src
OBJDIRS        =  build/x86_64-unknown-linux-gnu/support build/x86_64-unknown-linux-gnu/src

.PHONY: doc

all:
	@for DIR in $(OBJDIRS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    $(MAKE) ; status=$$? ; \
	    if [ $$status != 0 ] ; then \
	      echo "Exit status from make was $$status" ; exit $$status ; \
	    fi ; \
	    cd $$back ; \
	  done 
control:
	@for DIR in $(OBJDIRS) ;\
	  do LASTDIR=$$DIR ;\
	done ;\
	back=`pwd` ;\
	cd $$LASTDIR ;\
	$(MAKE) control_model  ;\
	cd $$back

one:
	@for DIR in $(OBJDIRS) ;\
	  do LASTDIR=$$DIR ;\
	done ;\
	back=`pwd` ;\
	cd $$LASTDIR ;\
	$(MAKE) $(name)  ;\
	cd $$back


install:
	@for DIR in $(OBJDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) install ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

clean:
	@for DIR in $(OBJDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	-rm -f ${exec_prefix}/bin/*
	-rm -f ${exec_prefix}/lib/*.a
	-rm -f ${exec_prefix}/module/*.mod  
	-rm -f ${exec_prefix}/src/*.o   
	-rm -rf html/[a-z]*

distclean:
	-rm -rf build
	-rm Makefile
	-rm build_command
	-rm config.log
	-rm config.status
	-rm config/config.h
	-rm config/mh-config-use
	-rm config/set-up.info
	-rm -rf doc/html
	-rm -rf doc/latex
	-rm -rf html/[a-z]*

doc:
	doxygen doc/resources/doxyfile_icon_html
	@echo 
	@echo "Start of HTML documentation: doc/html/index.html"
	@echo 

pdf: 
	doxygen doc/resources/doxyfile_icon_pdf

index:
	-rm -rf html/[a-z]*
	scripts/f2html_scripts/f2html.pl -f scripts/f2html_scripts/fgenrc -d html $(SRCDIRS)

checkstyle:
	scripts/codestyle_scripts/process_src -v
