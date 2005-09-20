##!
# \file         Makefile
# \author       Bill Hill
# \date         March 1999
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C) 2005 Medical research Council, UK.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief        Makefile for the Woolz type allocation library.
# \ingroup      Wlz
# \todo         -
# \bug          None known.
##

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			= libAlc.a
SOARCHIVE		= $(ARCHIVE:%.a=%.so)
DLLARCHIVE		= $(ARCHIVE:%.a=cyg%.dll)
DLLIMPORTLIB		= $(ARCHIVE:%.a=%.dll.a)
PROTOTYPES		= AlcProto.h
DLLENTRY		= dll_entry
DLLLIBDIRS		= $(DEPENDDIRS:%=-L%)
DLLLIBDEPS		= 
DLLLIBSEXT		=  

# Names of executables to be built (modify as required).
EXECUTABLES		=

# List of all 'C' source files (modify as required).
CSOURCES		= AlcAlloc.c \
			  AlcArray.c \
			  AlcBlockStack.c \
			  AlcCPQueue.c \
			  AlcDLPList.c \
			  AlcFreeStack.c \
			  AlcHashTable.c \
			  AlcKDTree.c \
			  AlcString.c \
			  AlcVector.c

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

# Basic flags for controlling compilation (modify as required).
DEFINES		= $(THREADS) \
		  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS	= $(INCDIRS:%=-I%) $(DEFINES)
CFLAGS		= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)
CXXFLAGS	= $(CXXDEBUG) $(CXXOPTIMISE) $(CXXANSI_CONFORMANCE)
 

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

# Make the tags file
tags:
			ctags $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)

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

AlcArray.o:		AlcArray.c AlcTemplates.h


# Targets to create executables, just a dependencies (modify as required).

# Target for Purified executables (modify as required).
purify:

