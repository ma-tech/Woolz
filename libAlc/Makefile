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
# Purpose:      Makefile for building the MRC HGU memory allocation
#		library.
# $Revision$
#########################################################################

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			= libAlc.a
SOARCHIVE		= $(ARCHIVE:%.a=%.so)
DLLARCHIVE		= $(ARCHIVE:lib%.a=%.dll)
DLLARCHIVEA		= $(ARCHIVE:%.a=%_dll.a)
DLLDEF			= $(ARCHIVE:%.a=%.def)
DLLBASE			= $(ARCHIVE:%.a=%.base)
DLLEXP			= $(ARCHIVE:%.a=%.exp)
PROTOTYPES		= AlcProto.h
DLLENTRY		= dll_entry
DLLLIBDEPS		= cygwin

# Names of executables to be built (modify as required).
EXECUTABLES		=

# List of all 'C' source files (modify as required).
CSOURCES		= AlcAlloc.c \
			  AlcArray.c \
			  AlcBlockStack.c \
			  AlcDLPList.c \
			  AlcFreeStack.c \
			  AlcHashTable.c \
			  AlcKDTree.c \
			  AlcString.c \
			  AlcVector.c

DLLCSOURCES		= init.c \
			  fixup.c

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		= Alc.h \
			  AlcType.h \
			  AlcProto.h

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		= AlcTemplates.h

# Source man documentation directory
SRC_MAN_DIR		= man

# Source html documentation directory
SRC_WWW_DIR		= html

# Destination html documentation directory
DST_WWW_DIR		= $(HGU_WWW_DIR)/src/Core/libAlc

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= .

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS)

# List of library search paths (modify as required).
LIBDIRS			= \
			  $(DEPENDDIRS)

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
DLLOBJECTS		= $(DLLCSOURCES:%.c=%.o)

# Basic flags for controlling compilation (modify as required).
DEFINES		= $(THREADS) \
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
			-$(RM) core $(OBJECTS) $(DLLOBJECTS) $(EXECUTABLES) \
			$(ARCHIVE) $(SOARCHIVE) \
			$(DLLARCHIVE) $(DLLBASE) \
			$(DLLARCHIVEA) $(DLLEXP) $(DLLDEF) \
			doxygen.tag *~
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

# Use doxygen to make the documentation (should not need modifying).
documentation:		doxygen.tag
			
doxygen.tag:		$(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)
			doxygen

# Make the tags file
tags:
			ctags $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)

# Targets to create archives, dependencies and rules (modify as required).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?
			$(RANLIB) $(ARCHIVE)

$(SOARCHIVE):		$(OBJECTS)
			$(CC) -G $(OBJECTS) -lC -o $(SOARCHIVE)

$(DLLARCHIVE):		$(OBJECTS) $(DLLOBJECTS) $(DLLDEF)
			ld -s --base-file $(DLLBASE) --dll -o $(DLLARCHIVE) \
			$(OBJECTS) $(DLLOBJECTS) $(DLLLIBDEPS:%=-l%) \
			-e _$(DLLENTRY)@12 $(DEPENDDIRS:%=-L%)
			dlltool --as=$(AS) --dllname $(DLLARCHIVE) \
			--def $(DLLDEF) --base-file $(DLLBASE) \
			--output-exp $(DLLEXP)
			ld -s --base-file $(DLLBASE) $(DLLEXP) \
			--dll -o $(DLLARCHIVE) \
			$(OBJECTS) $(DLLOBJECTS) $(DLLLIBDEPS:%=-l%) \
			-e _$(DLLENTRY)@12 $(DEPENDDIRS:%=-L%)
			dlltool --as=$(AS) --dllname $(DLLARCHIVE) \
			--def $(DLLDEF) --base-file $(DLLBASE) \
			--output-exp $(DLLEXP)
			ld $(DLLEXP) --dll -o $(DLLARCHIVE) $(OBJECTS) \
			$(DLLOBJECTS) $(DLLLIBDEPS:%=-l%) -e _$(DLLENTRY)@12 \
			$(DEPENDDIRS:%=-L%)
			rm -f $(DLLBASE) $(DLLEXP)
			dlltool --as=as --dllname $(DLLARCHIVE) \
			--def $(DLLDEF) --output-lib $(DLLARCHIVEA)

$(DLLDEF):		$(PROTOTYPES) $(OBJECTS) $(DLLOBJECTS)
			echo "EXPORTS" > $(DLLDEF)
			nm $(OBJECTS) $(DLLOBJECTS) | \
			grep '^........ [T] _' | \
			sed 's/[^_]*_//' >> $(DLLDEF)

# Targets to create executables, just a dependencies (modify as required).

# Target for Purified executables (modify as required).
purify:
