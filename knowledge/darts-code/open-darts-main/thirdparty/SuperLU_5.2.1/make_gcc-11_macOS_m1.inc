############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _mac_x

#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= .
SUPERLULIB   	= $(SuperLUroot)/libsuperlu_5.1.a

BLASLIB   	= $(SuperLUroot)/libblas.a

## ATLAS BLAS causes single-precision to fail
#BLASDEF = -DUSE_VENDOR_BLAS
#BLASLIB   	= -L/usr/lib/atlas -lblas

## This BLAS causes single-precision to fail the test in SuperLU
#BLASDEF 	= -DUSE_VENDOR_BLAS
#BLASLIB 	= -L/usr/lib -lblas

#TMGLIB       	= libtmglib.a
LIBS		= $(SUPERLULIB) $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = gcc-11
CFLAGS       = -O3
NOOPTS       =
FORTRAN	     = gfortran-11 #g77
FFLAGS       = -O2 -g -fopenmp
LOADER       = $(CC)
LOADOPTS     = -fopenmp

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
# This needs to be check if is needed or not
MATLAB	     = /Applications/MATLAB_R2022a.app
