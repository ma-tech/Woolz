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
ARCHIVE			= libJWlz.a
SOARCHIVE		= $(ARCHIVE:%.a=%.so)
DLLARCHIVE		= $(ARCHIVE:lib%.a=%.dll)
DLLARCHIVEA		= $(ARCHIVE:%.a=%_dll.a)
DLLDEF			= $(ARCHIVE:%.a=%.def)
DLLBASE			= $(ARCHIVE:%.a=%.base)
DLLEXP			= $(ARCHIVE:%.a=%.exp)
PROTOTYPES		= 
DLLENTRY		= dll_entry
DLLLIBDEPS		= Wlz_dll Alg_dll Alc_dll cygwin kernel32


# Names of executables to be built (modify as required).
EXECUTABLES		=

# List of all 'C' source files (modify as required).
CSOURCES		= $(shell ls Wlz*.c)
DLLCSOURCES		= init.c fixup.c

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
DLLOBJECTS		= $(DLLCSOURCES:%.c=%.o)

# Basic flags for controlling compilation (modify as required).
UNIXFLAGS		= -DUNIX
COPTIMISE		= -O2
ANSI_CONFORMANCE	=
DEFINES			= $(THREADS) \
			  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS		= $(INCDIRS:%=-I%) $(DEFINES)
CFLAGS			= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)

# Default traget which builds everything.
all:			archive

clean:
			$(RM) -f $(CSOURCES) *.h *.o *.so $(ARCHIVE) \
			$(DLLARCHIVE) $(DLLBASE) \
			$(DLLEXP) $(DLLDEF) *~

install_archive:	$(SOARCHIVE)
			$(INSTALL) -o $(OWNER) -m 644 $(SOARCHIVE) \
				   $(HGU_LIB_DIR)

archive:		$(ARCHIVE)

soarchive:		$(SOARCHIVE)

dll:			$(DLLARCHIVE)

$(ARCHIVE):		$(OBJECTS)
			$(AR) rs $(ARCHIVE) $?

$(SOARCHIVE):		$(OBJECTS) $(INCLUDES)
			$(CC) -G $(OBJECTS) $(LDFLAGS) $(LDLIBS) \
			      -o $(SOARCHIVE)

$(DLLARCHIVE):		$(OBJECTS) $(DLLOBJECTS) $(DLLDEF)
			dllwrap --output-def $(DLLDEF) --add-stdcall-alias \
			--driver-name $(CC) $(DEPENDDIRS:%=-L%) \
			-o $(DLLARCHIVE) $(OBJECTS) $(DLLOBJECTS) $(DLLLIBDEPS:%=-l%) \
			-W1,-e,_$(DLLENTRY)@12 -s

$(DLLDEF):		$(PROTOTYPES) $(OBJECTS) $(DLLOBJECTS)
			echo "EXPORTS" > $(DLLDEF)
			nm $(OBJECTS) $(DLLOBJECTS) | \
			grep '^........ [T] _' | \
			sed 's/[^_]*_//' >> $(DLLDEF)


