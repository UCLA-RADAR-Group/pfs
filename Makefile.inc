#
#	This is an include file which is included in all radar makefiles
#	Please edit this file
#	Please do not edit any of the other Makefiles
#
# You may have to indicate the correct path to sh for your system
SHELL   = /bin/sh
# 
# These are the C and FORTRAN compilers on your system 
CC      = gcc
F77     = g77
# 
# This is the library archiving command on your system
AR      = ar
#
# This is the name of your Fortran library
FLIB    = -lg2c
#
# This is the name and location of your X library
XLIB    = -L/usr/X11R6/lib -lX11
#
# This is where you would like the software to be installed 
DEST = $(HOME)
DESTDIR = $(DEST)/bin
DESTLIB = $(DEST)/lib
#
# Please do not change anything below this line
LIBDIR  = ../../lib
INCDIR  = ../../include
CFLAGS  = -I$(INCDIR) -O6 -DFFTW -DLARGEFILE -D_FILE_OFFSET_BITS=64
FFLAGS  = -I$(INCDIR)
LDFLAGS = -L../../lib -lm
BZFLAGS = -fno-second-underscore
