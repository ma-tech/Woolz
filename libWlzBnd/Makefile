#########################################################################
# Project:      Woolz
# Title:        Makefile
# Date:         March 1999
# Author:       Bill Hill
# Copyright:	2003 Medical Research Council, UK.
#		All rights reserved.
# Address:	MRC Human Genetics Unit,
#		Western General Hospital,
#		Edinburgh, EH4 2XU, UK.
# Purpose:      Makefile for the Woolz library.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			= libWlzBnd.a
SOARCHIVE		= $(ARCHIVE:%.a=%.so)
DLLARCHIVE		= $(ARCHIVE:%.a=cyg%.dll)
DLLIMPORTLIB		= $(ARCHIVE:%.a=%.dll.a)
PROTOTYPES		= WlzProto.h
DLLLIBDIRS		= $(DEPENDDIRS:%=-L%)
DLLLIBDEPS		= -lWlz -lAlg -lAlc
# DLLLIBSEXT		= -L. -lg -lm -lc
DLLLIBSEXT		= 

# Names of executables to be built (modify as required).
EXECUTABLES		= 

# List of all 'C' source files (modify as required).
CSOURCES		= \
			WlzBndGreyInvert.c

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		= \
			WlzBnd.h \
			WlzBndType.h \
			WlzBndProto.h

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		=

# Source man documentation directory
SRC_MAN_DIR		= man

# Source html documentation directory
SRC_WWW_DIR		= html

# Destination html documentation directory
DST_WWW_DIR		= $(HGU_WWW_DIR)/src/Core/libWlzBnd

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		=

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  . \
			  ../libAlc \
			  ../libAlg \
			  ../libWlz

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS) \
			  $(HGU_INC_DIR)

# List of library search paths (modify as required).
LIBDIRS			= \
			  $(DEPENDDIRS) \
			  $(SYSLIB) \
			  $(HGU_LIB_DIR)

# List of libraries to link (modify as required).
ifeq 		($(UNIXTYPE), SUNOS5)
EXTRA_LIBS		= gen socket nsl m
else
EXTRA_LIBS		= m
endif

MOTIFLIBS       	= Xm 
X11LIBS         	= Xt X11

LIBRARIES		=  WlzBnd Wlz Alg Alc \
			  $(EXTRA_LIBS)

# List of all sub directories that contain makefiles. The subsystems will be
# built in the order in which they are given here (modify as required).
SUBSYSTEMS		= 

# Name of distribution tar file (modify as required).

# List of files (should be top level directories) to be put used in creating
# a compressed tar file for distribution (modify as required).
TARSOURCES		= 

# List of all header files, ie the union of the private and public header
# files (should not need modifying).
INCLUDES_ALL		= $(sort $(INCLUDES_PUB) $(INCLUDES_PRV))

# List of all object files required for archive (should not need modifying).
OBJECTS			= $(CSOURCES:%.c=%.o)

# Basic flags for controlling compilation (modify as required).
DEFINES		= $(THREADS) \
		  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS	= $(INCDIRS:%=-I%) $(DEFINES)
CFLAGS		= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)
CXXFLAGS	=$(CXXDEBUG) $(CXXOPTIMISE) $(CXXANSI_CONFORMANCE)
ifeq ($(UNIXTYPE), SUNOS5)
ANSI_CONFORMANCE	=
endif

# Default traget which builds everything (should not need modifying).
all:			archive executables includes documentation
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# Target which installs everything after it has been built (should not need
# modifying).
install:		install_includes install_archive \
			install_executables install_documentation

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

# Target which installs the documentation in the appropriate directories,
# (should not need modifying).
install_documentation:	documentation
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_manpages)
			$(install_www)

# Target which installs the man pages in the appropriate directories,
# (should not need modifying).
install_manpages:	documentation
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			install_manpages

# Target which installs the html pages in the appropriate directories,
# (should not need modifying).
install_www:		documentation
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			install_www

# Target for compressed tar files suitable for distribution (should not need
# modifying).
tar:			$(TARSOURCES) clobber allsources
ifneq ($(strip $(TARSOURCES)),)
			$(TAR) $(TARFLAGS) $(TARFILE) $(TARSOURCES)
			$(COMPRESS) $(TARFILE)
endif

# Target for the archive (library), just a dependency (should not need
# modifying).
archive:		$(ARCHIVE)

soarchive:		$(SOARCHIVE)

dll:			$(DLLARCHIVE)

# Target for executables, just a dependency (should not need modifying).
executables:		$(EXECUTABLES)

# Target for public header files, just a dependency (should not need modifying).
includes:		$(INCLUDES_PUB)

# Target for manual pages, just a dependency (should not need modifying).
manpages:		$(MANPAGES)

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

# Use doxygen to make the documentation (should not need modifying).
documentation:		doxygen.tag

doxygen.tag:		$(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)
			doxygen

# making the tag file
tags:			 $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)
			 $(CTAGS) $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)

allsources:		sources
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

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

# Target for Purified executables (modify as required).
purify:

$(OBJECTS):		%.o: %.c $(INCLUDES_ALL) 
			$(CC) -c -o $@ $(CPPFLAGS) $(CFLAGS) $<

