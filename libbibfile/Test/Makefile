#########################################################################
# Project:      Mouse Atlas
# Title:        Makefile
# Date:		March 1999
# Author:	Bill Hill
# Copyright:	1999 Medical Research Council, UK.
#		All rights reserved.
#		MRC Human Genetics Unit,
#		Western General Hospital,
#		Edinburgh, EH4 2XU, UK.
# Purpose:      Makefile for the bibfile test executable.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			=

# Names of executables to be built (modify as required).
EXECUTABLES		= bibFileTest

# List of all 'C' source files (modify as required).
CSOURCES		= bibFileTest.c

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		=

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		= bibFile.h

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		=

# List of all header file directories to be searched (modify as required).
INCDIRS			= ../  \
			  $(HGUHOME)/include

# List of library search paths (modify as required).
LIBDIRS			= .. \
			  $(SYSLIB)

# List of libraries to link (modify as required).
LIBRARIES		= bibfile \
			  Alc \
			  gen


# Linker library flags (should not need modifying).
LDLIBS			= $(LIBRARIES:%=-l%)

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
DEFINES			= -D$(UNIXTYPE)
CPPFLAGS		= $(INCDIRS:%=-I%) $(DEFINES)
#CDEBUG			=
CDEBUG			= -g -DDEBUG
COPTIMISE		= -KPIC -xchip=ultra
#COPTIMISE		= 
CFLAGS			= $(CDEBUG) $(COPTIMISE)

# Default traget which builds everything (should not need modifying).
all:			archive executables includes manpages
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# Target which installs everything after it has been built (should not need
# modifying).
install:		install_includes install_manpages install_archive \
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

# Target which installs the manual pages in the appropriate directories,
# watch out for the tabs in the sed regexp's (should not need modifying).
install_manpages:	$(MANPAGES)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_manpages)

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

# Target for manual pages, just a dependency (should not need modifying).
manpages:		$(MANPAGES)

# Target for lint (modify as required).
lint:			$(CSOURCES) $(INCLUDES_PRV)
			$(LINT) $(LINTFLAGS) $(CPPFLAGS) $(CSOURCES)

# Target to clean up after making (should not need modifying).
clean:			
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

# Target to create archive (library) (should not need modifying).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?

# Targets to create executables, just a dependencies (modify as required).
bibFileTest:		bibFileTest.o

# Target for Purified executables, just a dependencies (modify as required).
purify:			bibFileTest.pure
bibFileTest.pure:	bibFileTest.o
