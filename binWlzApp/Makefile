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
# Purpose:      Makefile for small non-window system applications
#		which use Woolz and other Core libraries.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			=

# Names of executables to be built (modify as required).
EXECUTABLES		= \
			Wlz3DGetSection \
			Wlz3DViewTransformObj \
			Wlz3DViewTransformBitmap \
			WlzExtractTransform \
			WlzProductTransform \
			WlzFixedPlaneAlign \
			WlzPatchObjRegister \
			WlzRawToWlz \
			WlzSizeSelect \
			WlzSplitImage \
			WlzObjCompareSpecial_01 \
			WlzVtxToObj \
			WlzFixDomains \
			WlzMatchICPPlane \
			WlzTiffStackToShade \
			WlzTiffStackToPatches

SCRIPTS			= \
			WlzObjCompareSpecial_01.csh

# List of all 'C' source files (modify as required).
CSOURCES		= \
			Wlz3DGetSection.c \
			Wlz3DViewTransformObj.c \
			Wlz3DViewTransformBitmap.c \
			WlzExtractTransform.c \
			WlzProductTransform.c \
			WlzFixedPlaneAlign.c \
			WlzPatchObjRegister.c \
			WlzRawToWlz.c \
			WlzSizeSelect.c \
			WlzSplitImage.c \
			WlzObjCompareSpecial_01.c \
			WlzVtxToObj.c \
			WlzFixDomains.c \
			WlzMatchICPPlane.c \
			WlzTiffStackToShade.c \
			WlzTiffStackToPatches.c

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
			  ../libAlc \
			  ../libAlg \
			  ../libWlz \
			  ../libbibfile \
			  ../libWlzExtFF \
			  ../libhguDlpList \
			  ../libReconstruct \
			  ../../External/Tiff/tiff-v3.5.7/libtiff \
			  ../../External/Jpeg/jpeg-6b 

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS) \
			  $(HGU_INC_DIR) \
			  /opt/local/include

# List of library search paths (modify as required).
LIBDIRS			= \
			  $(DEPENDDIRS) \
			  $(HGU_LIB_DIR) \
			  $(SYSLIB)

# List of libraries to link (modify as required).
ifeq 		($(UNIXTYPE), SUNOS5)
EXTRA_LIBS		= jpeg tiff gen socket nsl m
else
EXTRA_LIBS		= jpeg tiff m
endif

MOTIFLIBS       	= Xm 
X11LIBS         	= Xt X11

LIBRARIES		=  Reconstruct WlzExtFF hguDlpList bibfile \
			   Wlz Alg Alc \
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
# DEFINES			= -DUSE_OSTORE -D__EXTENSIONS__
DEFINES			= 
CPPFLAGS		= $(INCDIRS:%=-I%) $(DEFINES) $(UNIXFLAGS)
CFLAGS			= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)

# Default traget which builds everything (should not need modifying).
all:			archive executables includes manpages
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# Target which installs everything after it has been built (should not need
# modifying).
install:		install_includes install_manpages install_archive \
			install_executables install_release

# Target which installs an archive (should not need modifying).
install_archive:	$(ARCHIVE)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_archive)

# Target which installs executables (should not need modifying).
install_executables:	$(EXECUTABLES) $(SCRIPTS)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_executables)

# Target which installs release executables (should not need modifying).
install_release:	$(EXECUTABLES) $(SCRIPTS)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_release)


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
tar:			$(TARSOURCES) clobber
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

# Targets to create archives, dependencies and rules (modify as required).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?

# Targets to create executables, just a dependencies (modify as required).

# Target for Purified executables (modify as required).
purify:

$(OBJECTS):		%.o: %.c $(INCLUDES_ALL) 
ifeq ($(UNIXTYPE), SUNOS5)
			$(CC) -v  -c -o $@ $(CPPFLAGS) $(CFLAGS) $<
else
			$(CC) -c -o $@ $(CPPFLAGS) $(CFLAGS) $<
endif
