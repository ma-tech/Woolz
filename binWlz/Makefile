#########################################################################
# Project:      Woolz
# Title:        Makefile
# Date:         September 1994
# Author:       Bill Hill
# Copyright:	1999 Medical Research Council, UK.
#		All rights reserved.
# Address:	MRC Human Genetics Unit,
#		Western General Hospital,
#		Edinburgh, EH4 2XU, UK.
# Purpose:      Makefile for the Woolz filters.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
# 07-02-2k bill	Add WlzHistogramFitPeaks.
# 02-02-2k bill	Add WlzHistogramFindPeaks and WlzHistogramGauss.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			=

# Names of executables to be built (modify as required).
EXECUTABLES		= \
			Wlz3DGetSection \
			Wlz3DViewTransformObj \
			WlzAffineTransformLSq \
			WlzAffineTransformObj \
			WlzArea \
			WlzBasisFnTransformObj \
			WlzBoundaryToObj \
			WlzBoundingBox \
			WlzCannyDeriche \
			WlzCentreOfMass \
			WlzClipObjToBox \
			WlzCompThresh \
			WlzCompound \
			WlzConvertPix \
			WlzConvexPoly \
			WlzCutObjToBox \
			WlzDilation \
			WlzDiffDomain \
			WlzDomain \
			WlzDomainFill \
			WlzErosion \
			WlzExplode \
			WlzFacts \
			WlzGauss \
			WlzGaussNoise \
			WlzGreyCrossing \
			WlzGreyGradient \
			WlzGreyInvert \
			WlzGreyMask \
			WlzGreyNormalise \
			WlzGreyRange \
			WlzGreySetRange \
			WlzGreySetValue \
			WlzGreyStats \
			WlzGreyTemplate \
			WlzHistogramData \
			WlzHistogramEqualiseObj \
			WlzHistogramFindPeaks \
			WlzHistogramFitPeaks \
			WlzHistogramGauss \
			WlzHistogramMatchObj \
			WlzHistogramObj \
			WlzHistogramRebin \
			WlzHistogramSmooth \
			WlzImageArithmetic \
			WlzIntersect \
			WlzLabel \
			WlzLaplacian \
			WlzObjToBoundary \
			WlzPolarSample \
			WlzPolygonToObj \
			WlzPrinicipalAngle \
			WlzRsvFilterObj \
			WlzSampleObj \
			WlzSelect1InN \
			WlzSetBackground \
			WlzSetVoxelSize \
			WlzShiftObj \
			WlzSkeleton \
			WlzSobel \
			WlzTransposeObj \
			WlzThreshold \
			WlzUnion \
			WlzVerifyObj \
			WlzVolume \
			WlzWindow

#			WlzAutoTrackUpDown


# List of all 'C' source files (modify as required).
CSOURCES		= \
			MADomainTest.c \
			Wlz3DGetSection.c \
			WlzAffineTransformLSq.c \
			WlzAffineTransformObj.c \
			WlzArea.c \
			WlzAutoTrackUpDown.c \
			WlzBoundaryToObj.c \
			WlzBoundingBox.c \
			WlzCannyDeriche.c \
			WlzCentreOfMass.c \
			WlzClipObjToBox.c \
			WlzCompound.c \
			WlzCompThresh.c \
			WlzConvertPix.c \
			WlzConvexPoly.c \
			WlzCutObjToBox.c \
			WlzDilation.c \
			WlzDiffDomain.c \
			WlzDomain.c \
			WlzDomainFill.c \
			AlgDPSearch.c \
			WlzErosion.c \
			WlzExplode.c \
			WlzExplode3D.c \
			WlzFacts.c \
			WlzGauss.c \
			WlzGaussNoise.c \
			WlzGreyCrossing.c \
			WlzGreyGradient.c \
			WlzGreyInvert.c \
			WlzGreyMask.c \
			WlzGreyNormalise.c \
			WlzGreyRange.c \
			WlzGreySetRange.c \
			WlzGreySetValue.c \
			WlzGreyStats.c \
			WlzGreyTemplate.c \
			WlzHistogramData.c \
			WlzHistogramEqualiseObj.c \
			WlzHistogramFindPeaks.c \
			WlzHistogramFitPeaks.c \
			WlzHistogramGauss.c \
			WlzHistogramMatchObj.c \
			WlzHistogramObj.c \
			WlzHistogramRebin.c \
			WlzHistogramSmooth.c \
			WlzImageArithmetic.c \
			WlzIntersect.c \
			WlzLabel.c \
			WlzLaplacian.c \
			WlzObjToBoundary.c \
			WlzPolarSample.c \
			WlzPolygonToObj.c \
			WlzPrinicipalAngle.c \
			WlzRsvFilterObj.c \
			WlzSampleObj.c \
			WlzSelect1InN.c \
			WlzSetBackground.c \
			WlzSetVoxelSize.c \
			WlzShiftObj.c \
			WlzSkeleton.c \
			WlzSobel.c \
			WlzTransposeObj.c \
			WlzThreshold.c \
			WlzUnion.c \
			WlzVerifyObj.c \
			WlzVolume.c \
			WlzWindow.c

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		= 

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		=

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		= \
			  WlzAffineTransformLSq.1 \
			  WlzAffineTransformObj.1 \
			  WlzBoundingBox.1 \
			  WlzCentreOfMass.1 \
			  WlzClipObjToBox.1 \
			  WlzCutObjToBox.1 \
			  WlzErosion.1 \
			  WlzExplode.1 \
			  WlzFacts.1 \
			  WlzGreyCrossing.1 \
			  WlzGreyStats.1 \
			  WlzHistogramData.1 \
			  WlzHistogramEqualiseObj.1 \
			  WlzHistogramMatchObj.1 \
			  WlzHistogramObj.1 \
			  WlzHistogramRebin.1 \
			  WlzHistogramSmooth.1 \
			  WlzLaplacian.1 \
			  WlzMorphological.1 \
			  WlzPolarSample.1 \
			  WlzPrinicipalAngle.1 \
			  WlzSampleObj.1 \
			  WlzShiftObj.1 \
			  WlzSobel.1 \
			  WlzWindow.1

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  . \
			  ../libAlc \
			  ../libAlg \
			  ../libbibfile \
			  ../libWlz

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS) \
			  $(HGU_INC_DIR) \
                          $(OS_ROOTDIR)/include 

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

LIBRARIES		= Wlz bibfile Alg Alc \
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
#CDEBUG			= -g 
#COPTIMISE		= -xchip=ultra -KPIC
#CXXDEBUG		= -g 
#CXXOPTIMISE		= -xchip=ultra -KPIC

# un-comment these lines to make a 32-bit library for the IRIX 6
ifeq ($(UNIXTYPE), IRIX5)
ANSI_CONFORMANCE	= -xansi -fullwarn
CC			= cc -32
COPTIMISE		= -mips2 -O1
SYSLIB			= /usr/lib
UNIXFLAGS		= -DUNIX  -DSYSV -DSVR4 
endif	

ifeq 		($(UNIXTYPE), IRIX646)
CDEBUG			= 
COPTIMISE		= -O1 -OPT:roundoff=3
endif

DEFINES			= -DUSE_OSTORE
#DEFINES			= 
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

# Targets to create archives, dependencies and rules (modify as required).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?

# Targets to create executables, just a dependencies (modify as required).
WlzAutoTrackUpDown:	WlzAutoTrackUpDown.o AlgDPSearch.o
WlzObjToBoundary:	WlzObjToBoundary.o

# Target for Purified executables (modify as required).
purify:

$(OBJECTS):		%.o: %.c $(INCLUDES_ALL) 
ifeq ($(UNIXTYPE), SUNOS5)
			$(CC) -v  -c -o $@ $(CPPFLAGS) $(CFLAGS) $<
else
			$(CC) -c -o $@ $(CPPFLAGS) $(CFLAGS) $<
endif
