#########################################################################
# Project:      Mouse Atlas
# Title:        Makefile
# Date:         March 1999
# Author:       Bill Hill
# Copyright:	1999 Medical Research Council, UK.
#		All rights reserved.
# Address:	MRC Human Genetics Unit,
#		Western General Hospital,
#		Edinburgh, EH4 2XU, UK.
# Purpose:      Makefile for building the MRC Human Genetics Unit
#		numerical algorithm library.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
# 07-06-01 bill	Add AlgMatrixMath.c
# 04-06-01 bill Add AlgMatrixRSTDiag.c, AlgMatrixTDiagQLI.c and
#		AlgMatrixRSEigen.c
# 01-02-01 bill	Add AlgCrossCorr.c
# 07-11-00 bill Add AlgShuffle.c.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			= libAlg.a
SOARCHIVE		= $(ARCHIVE:%.a=%.so)
DLLARCHIVE		= $(ARCHIVE:lib%.a=cyg%.dll)
DLLIMPORTLIB		= $(ARCHIVE:%.a=%.dll.a)
PROTOTYPES		= AlgProto.h
DLLLIBDIRS		= $(DEPENDDIRS:%=-L%)
DLLLIBDEPS		= -lAlc
DLLLIBSEXT		= 

# Names of executables to be built (modify as required).
EXECUTABLES		=

# List of all 'C' source files (modify as required).
CSOURCES		= \
			  AlgBits.c \
			  AlgComplexUtils.c \
			  AlgConvolve.c \
			  AlgCrossCorr.c \
			  AlgDebug.c \
			  AlgDPSearch.c \
			  AlgFourier.c \
			  AlgGamma.c \
			  AlgHeapSort.c \
			  AlgLinearFit.c \
			  AlgMatrixCG.c \
			  AlgMatrixGauss.c \
			  AlgMatrixLU.c \
			  AlgMatrixMath.c \
			  AlgMatrixRSEigen.c \
			  AlgMatrixRSTDiag.c \
			  AlgMatrixSV.c \
			  AlgMatrixTDiagQLI.c \
			  AlgMixture.c \
			  AlgPolyLSQ.c \
			  AlgQSort.c \
			  AlgRand.c \
			  AlgRange.c \
			  AlgRank.c \
			  AlgShuffle.c \
			  AlgVectorMath.c
			  
# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		= Alg.h \
			  AlgType.h \
			  AlgProto.h

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		= Alg.h \
			  AlgType.h \
			  AlgProto.h

# Source man documentation directory
SRC_MAN_DIR		= man

# Source html documentation directory
SRC_WWW_DIR		= html

# Destination html documentation directory
DST_WWW_DIR		= $(HGU_WWW_DIR)/src/Core/libAlg

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		= 

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= . \
			  ../libAlc

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS)  \
			  $(HGUHOME)/include

# List of library search paths (modify as required).
LIBDIRS			=

# List of all sub directories that contain makefiles (modify as required).
SUBSYSTEMS		=

# List of files (should be top level directories) to be put used in creating
# a compressed tar file for distribution.
TARSOURCES		=

# List of all header files, ie the union of the private and public header
# files (should not need modifying).
INCLUDES_ALL		= $(sort $(INCLUDES_PUB) $(INCLUDES_PRV))

# List of all object files required for archive (should not need modifying).
OBJECTS			= $(CSOURCES:%.c=%.o)

# Basic flags for controlling compilation (modify as required).
DEFINES		= $(THREADS) -D__EXTENSIONS__\
		  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS	= $(INCDIRS:%=-I%) $(DEFINES)
CFLAGS		= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)
CXXFLAGS	=$(CXXDEBUG) $(CXXOPTIMISE) $(CXXANSI_CONFORMANCE)
 

# Default traget which builds everything (should not need modifying).
all:			archive executables includes documentation
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# Target which installs everything after it has been built (should not need
# modifying).
install:		install_includes install_archive install_executables \
			install_documentation

# Target which installs an archive (should not need modifying).
install_archive:	$(ARCHIVE)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_archive)

# Target which installs executables (should not need modifying).
install_executables:	$(EXECUTABLES)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_executables)

# Target which installs the public 'C' header files (should not need 
# modifying).
install_includes:	$(INCLUDES_PUB)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_includes)

# Target for the archive (library), just a dependency (should not need
# modifying).
archive:		$(ARCHIVE)

soarchive:		$(SOARCHIVE)

dll:			$(DLLARCHIVE)

# Target for executables, just a dependency (should not need modifying).
executables:		$(EXECUTABLES)

# Target for public header files, just a dependency (should not need modifying).
includes:		$(INCLUDES_PUB)

# Target for lint (modify as required).
lint:			$(CSOURCES) $(INCLUDES_PRV)
			$(LINT) $(LINTFLAGS) $(CPPFLAGS) $(CSOURCES)

# Target to clean up after making (should not need modifying).
clean:
			-$(RM) core $(OBJECTS) $(EXECUTABLES) \
			$(ARCHIVE) $(SOARCHIVE) $(DLLARCHIVE) \
			*.dll.a doxygen.tag *~
			-$(RMDIR) man html

# Target to clean up everything that might not be wanted (should not need
# modifying).
clobber:		clean
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			-$(RM) *.o *.a .pure *.pure *.bak *~

allsources:		sources
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# making the tag file
tags:
	$(CTAGS)  $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)

# Targets to create archives, dependencies and rules (modify as required).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?
			$(RANLIB) $(ARCHIVE)

$(SOARCHIVE):		$(OBJECTS)
			$(CC) -G $(OBJECTS) -lC -o $(SOARCHIVE)


$(DLLARCHIVE):		$(OBJECTS)
			$(CC) -shared -o $(DLLARCHIVE) \
			-Wl,--out-implib=$(DLLIMPORTLIB) \
			-Wl,--export-all-symbols \
			-Wl,--enable-auto-import \
			-Wl,--whole-archive $(OBJECTS) $(DLLLIBDIRS) $(DLLLIBDEPS) \
			-Wl,--no-whole-archive $(DLLLIBSEXT)

# Targets to create executables, just a dependencies (modify as required).
$(OBJECTS):		%.o: %.c $(INCLUDES_ALL) 

# Target for Purified executables (modify as required).
purify:



