#########################################################################
# Project:      Woolz
# Title:        Makefile
# Date:         August 1999
# Author:       Bill Hill
# Copyright:	1999 Medical Research Council, UK.
#		All rights reserved.
# Address:	MRC Human Genetics Unit,
#		Western General Hospital,
#		Edinburgh, EH4 2XU, UK.
# Purpose:      Makefile for the Java/Woolz binding.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../../../../../../../Makefile.conf

# Name of archive (library) to be built (modify as required).
SOARCHIVE		= libJWlz.so

# Names of executables to be built (modify as required).
EXECUTABLES		=

# List of all 'C' source files (modify as required).
CSOURCES		= $(shell ls *.c)

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		=

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		= $(shell ls *.h)

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		=

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  . \
			  ../../../../../../../libAlc \
			  ../../../../../../../libAlg \
			  ../../../../../../../libWlz

# List of all header file directories to be searched (modify as required).
INCDIRS                 = \
			  $(DEPENDDIRS) \
			  $(HGU_INC_DIR) \
			  $(JAVA_INC_DIR)

# List of library search paths (modify as required).
LIBDIRS			= \
			  $(DEPENDDIRS) \
			  $(SYSLIB) \
			  $(HGU_LIB_DIR)

# List of libraries to link (modify as required).
EXTRA_LIBS		= gen socket nsl m
LIBRARIES		= Wlz Alg Alc \
			  $(EXTRA_LIBS)

# List of all sub directories that contain makefiles. The subsystems will be
# built in the order in which they are given here (modify as required).
SUBSYSTEMS		=

# List of all object files required for archive (should not need modifying).
OBJECTS			= $(CSOURCES:%.c=%.o)

# Basic flags for controlling compilation (modify as required).
DEFINES			= $(THREADS) \
			  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS		= $(INCDIRS:%=-I%) $(DEFINES)
CFLAGS			= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)

# Default traget which builds everything.
all:			archive

clean:
			$(RM) *.o *.so *.c *.h

install_archive:	$(SOARCHIVE)
			$(INSTALL) -o $(OWNER) -m 644 $(SOARCHIVE) \
				   $(HGU_LIB_DIR)
			

archive:		$(SOARCHIVE)

$(SOARCHIVE):		$(OBJECTS) $(INCLUDES)
			$(CC) -G $(OBJECTS) $(LDFLAGS) $(LDLIBS) \
			      -o $(SOARCHIVE)
