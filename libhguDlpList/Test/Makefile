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
# Purpose:      Makefile which builds test executable for the HGU
#		doubly linked list of pointers library.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			=

# Names of executables to be built (modify as required).
EXECUTABLES		= HGUDlpListTest

# List of all 'C' source files (modify as required).
CSOURCES		= HGUDlpListTest.c

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		= 

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		= 

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		=

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  . \
			  ..

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS) \
			  $(HGUHOME)/include

# List of library search paths (modify as required).
LIBDIRS			= \
			  $(DEPENDDIRS) \
			  $(HGUHOME)/lib

# List of libraries to link (modify as required).
LIBRARIES		= hguDlpList

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
DEFINES                 = 
CPPFLAGS                = $(INCDIRS:%=-I%) $(DEFINES) $(UNIXFLAGS)
CFLAGS                  = $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)

# Default traget which builds everything (should not need modifying).
all:			archive executables includes manpages subsystems

# Target which installs everything after it has been built (should not need
# modifying).
install:		install_archive install_executables install_includes \
			install_manpages

# Target which installs an archive (should not need modifying).
install_archive:	$(ARCHIVE)
ifneq ($(strip $(ARCHIVE)),)
			$(INSTALL) -m 644 -o root $(ARCHIVE) $(HGU_LIB_DIR)
endif

# Target which installs executables (should not need modifying).
install_executables:	$(EXECUTABLES)
ifneq ($(strip $(EXECUTABLES)),)
			$(INSTALL) -m 644 -o root $(EXECUTABLES) $(HGU_BIN_DIR)
endif

# Target which installs the public 'C' header files (should not need 
# modifying).
install_includes:
ifneq ($(strip $(INCLUDES_PUB)),)
			$(INSTALL) -m 644 -o root $(INCLUDES_PUB) $(HGU_INC_DIR)
endif

# Target which installs the manual pages in the appropriate directories,
# watch out for the tabs in the sed regexp's (should not need modifying).
install_manpages:	$(MANPAGES)
ifneq ($(strip $(MANPAGES)),)
			for N in 1 2 3 4 5 6 7 8 9 ; \
			do \
			  M=`echo $(MANPAGES) | \
			     sed -e "s/[^. 	][^. 	]*\\.[^$$N]/ /g"` ; \
			  if [ -n "$$M" ] ; \
			  then \
			    for P in $$M ; \
			    do \
			      $(INSTALL) -m 644 -o root $$P \
					 $(HGU_MAN_DIR)/man.$$N ; \
			    done \
			  fi \
			done
endif

# Target for compressed tar files suitable for distribution (should not need
# modifying).
tar:			$(TARSOURCES)
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

# Target for making targets in sub directories (should not need modifying).
subsystems:		
ifneq ($(strip $(SUBSYSTEMS)),)
			for SUBSYSTEM in $(SUBSYSTEMS) ; \
			do  \
			  $(MAKE) -C $$SUBSYSTEM ; \
			done
endif

# Target to clean up after making (should not need modifying).
clean:
			-$(RM) $(OBJECTS) $(EXECUTABLES) $(ARCHIVE)

# Target to clean up everything that might not be wanted (should not need
# modifying).
clobber:		clean
			-$(RM) *.o *.a .pure *.pure *.bak *~
ifneq ($(strip $(SUBSYSTEMS)),)
			for SUBSYSTEM in $(SUBSYSTEMS) ; \
			do  \
			  $(MAKE) -C $$SUBSYSTEM clobber ; \
			done
endif

# Target to create archive (library) (should not need modifying).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?

# Targets to create executables, just a dependencies (modify as required).
HGUDlpListTest:		HGUDlpListTest.o

# Target for Purified executables (modify as required).
purify:			HGUDlpListTest.pure
HGUDlpListTest.pure:	HGUDlpListTest.o
