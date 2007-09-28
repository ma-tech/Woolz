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
# \brief        Makefile for the bibtex parsing library with functions
#		for parsing bibtex style text files as used by the Mouse
#		Atlas Project.
# \ingroup      Alg
# \todo         -
# \bug          None known.
##

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			= libbibfile.a

# Names of executables to be built (modify as required).
EXECUTABLES		=

# List of all 'C' source files (modify as required).
CSOURCES		= bibFileIO.c \
			  bibFileParse.c \
			  bibFileAlloc.c

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		= bibFile.h

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		= bibFile.h

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		= bibfile.3 \
			  bibfile.5

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  . \
			  ../libAlc

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS)/  \
			  $(HGU_INC_DIR)/include

# List of library search paths (modify as required).
LIBDIR			=

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
CXXFLAGS	=$(CXXDEBUG) $(CXXOPTIMISE) $(CXXANSI_CONFORMANCE)
 

# Default traget which builds everything (should not need modifying).
all:			archive executables includes tag
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# Target which installs everything after it has been built (should not need
# modifying).
install:		install_includes install_archive \
			install_executables

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

# Target for executables, just a dependency (should not need modifying).
executables:		$(EXECUTABLES)

# Target for public header files, just a dependency (should not need modifying).
includes:		$(INCLUDES_PUB)

# Target for tags files (modify as required).
tag:			tags
 
tags:			$(CSOURCES) $(INCLUDES_ALL)
			$(CTAGS) $(CSOURCES) $(INCLUDES_ALL)
 
# Target for lint (modify as required).
lint:			$(CSOURCES) $(INCLUDES_PRV)
			$(LINT) $(LINTFLAGS) $(CPPFLAGS) $(CSOURCES)

# Target to clean up after making (should not need modifying).
clean:
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			-$(RM) core $(OBJECTS) $(EXECUTABLES) $(ARCHIVE)

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

# Targets to create archives, dependencies and rules (modify as required).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?
			$(RANLIB) $(ARCHIVE)

# Targets to create executables, just a dependencies (modify as required).

# Target for Purified executables (modify as required).
purify:
