#########################################################################
# Project:      Woolz
# Title:        Makefile
# Date:         March 1999
# Author:       Bill Hill
# Copyright:	1999 Medical Research Council, UK.
#		All rights reserved.
# Address:	MRC Human Genetics Unit,
#		Western General Hospital,
#		Edinburgh, EH4 2XU, UK.
# Purpose:      Makefile for the Woolz library.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
# 29-11-99 bill	Add WlzContour.c
#########################################################################

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			= libWlz.a
SOARCHIVE		= $(ARCHIVE:%.a=%.so)

# Names of executables to be built (modify as required).
EXECUTABLES		= 

# List of all 'C' source files (modify as required).
CSOURCES		= \
			Wlz3DSection.c \
			Wlz3DSectionSegmentObject.c \
			Wlz3DViewStructUtils.c \
			Wlz3DViewTransformObj.c \
			WlzAffineTransform.c \
			WlzAffineTransformLSq.c \
			WlzArea.c \
			WlzArray.c \
                        WlzAssign.c \
			WlzBackground.c \
			WlzBasisFnTransform.c \
			WlzBoundingBox.c \
			WlzBoundToObj.c \
			WlzCannyDeriche.c \
			WlzCentreOfMass.c \
			WlzClipObjToBox.c \
			WlzCompThresh.c \
			WlzContour.c \
			WlzConvertPix.c \
			WlzConvexHull.c \
			WlzConvolve.c \
			WlzCopy.c \
			WlzCutObjToBox.c \
			WlzDebug.c \
			WlzDiffDomain.c \
			WlzDiffDomain3d.c \
			WlzDilation.c \
			WlzDomainFill.c \
			WlzDomainUtils.c \
			WlzErosion.c \
			WlzExplode3D.c \
			WlzErosion4.c \
			WlzError.c \
			WlzFacts.c \
			WlzFillBlankPlanes.c \
			WlzFreeSpace.c \
			WlzGauss.c \
			WlzGaussNoise.c \
			WlzGeometry.c \
			WlzGreyCrossing.c \
			WlzGreyDitherObj.c \
			WlzGreyGradient.c \
			WlzGreyInvertMinMax.c \
			WlzGreyMask.c \
			WlzGreyModGradient.c \
			WlzGreyNormalise.c \
			WlzGreyRange.c \
			WlzGreyScan.c \
			WlzGreySetRange.c \
			WlzGreySetValue.c \
			WlzGreyStats.c \
			WlzGreyTemplate.c \
			WlzGreyValue.c \
			WlzHistogram.c \
			WlzHyThreshold.c \
			WlzImageArithmetic.c \
			WlzInsideDomain.c \
			WlzIntersect2.c \
			WlzIntersect3d.c \
			WlzIntersectN.c \
			WlzIntervalCount.c \
			WlzIntervalDomScan.c \
			WlzIntRescaleObj.c \
			WlzLabel.c \
			WlzLabel3d.c \
			WlzLaplacian.c \
			WlzLineArea.c \
			WlzMakeAffineTransform.c \
			WlzMakeCompound.c \
			WlzMakeIntervalValues.c \
			WlzMakeProperties.c \
			WlzMakeStructs.c \
			WlzMeshTransform.c \
			WlzMeshUtils.c \
			WlzNMSuppress.c \
			WlzObjToBoundary.c \
			WlzPolarSample.c  \
			WlzPolyToObj.c \
			WlzPrinicipalAngle.c \
			WlzReadObj.c \
			WlzRsvFilter.c \
			WlzSampleObj.c \
			WlzScalarArithmeticOp.c \
			WlzScalarBinaryOp.c \
			WlzSepTrans.c \
			WlzSeqPar.c \
			WlzShift.c \
			WlzSkeleton.c \
			WlzSobel.c \
			WlzStdStructElements.c \
			WlzStringTypes.c \
			WlzStringUtils.c \
			WlzStructDilation.c \
			WlzStructErosion.c \
			WlzThreshold.c \
			WlzTransposeObj.c \
			WlzUnion2.c \
			WlzUnion3d.c \
			WlzUnionN.c \
			WlzValueTableUtils.c \
			WlzValueUtils.c \
			WlzVerifyObj.c \
			WlzVolume.c \
			WlzWindow.c \
			WlzWriteObj.c

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		= \
			Wlz.h \
			WlzType.h \
			WlzMacro.h \
			WlzError.h \
			WlzDebug.h \
			WlzProto.h

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		=

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		= \
			Wlz3DSectionIncrementDistance.3 \
			Wlz3DViewUtils.3 \
			WlzAffineTransform.3 \
			WlzAffineTransformLSq.3 \
			WlzArithmeticOps.3 \
			WlzArray.3 \
			WlzAssign.3 \
			WlzBackground.3 \
			WlzBinaryOps.3 \
			WlzBoundingBox.3 \
			WlzCannyDeriche.3 \
			WlzCentreOfMass.3 \
			WlzClipObjToBox.3 \
			WlzConvolve.3 \
			WlzCopy.3 \
			WlzCutObjToBox.3 \
			WlzDebug.3 \
			WlzDomainFill.3 \
			WlzDomainUtils.3 \
			WlzExplode3D.3 \
			WlzFileIO.3 \
			WlzFree.3 \
			WlzGetSectionFromObject.3 \
			WlzGreyCrossing.3 \
			WlzGreyFunctions.3 \
			WlzGreyGradient.3 \
			WlzGreyRange.3 \
			WlzGreyScan.3 \
			WlzGreyStats.3 \
			WlzGreyTableUtils.3 \
			WlzGreyValue.3 \
			WlzHistogram.3 \
			WlzHyThreshold.3 \
			WlzInsideDomain.3 \
			WlzIntervalScan.3 \
			WlzLaplacian.3 \
			WlzMakeSpace.3 \
			WlzMorphologicalOps.3 \
			WlzMorphometry.3 \
			WlzNMSuppress.3 \
			WlzObjectFacts.3 \
			WlzPolarSample.3 \
			WlzPrinicipalAngle.3 \
			WlzRsvFilter.3 \
			WlzSampleObj.3 \
			WlzSegment.3 \
			WlzSeqPar.3 \
			WlzShift.3 \
			WlzSobel.3 \
			WlzTranspose.3 \
			WlzValueUtils.3 \
			WlzVerify.3 \
			WlzWindow.3

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  . \
			  ../libAlc \
			  ../libAlg

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

LIBRARIES		=  Wlz Alg Alc \
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
#CDEBUG		= -g -DDEBUG 
#COPTIMISE	= -xchip=ultra -KPIC
#CXXDEBUG	= -g -DDEBUG
#CXXOPTIMISE	= -xchip=ultra -KPIC
DEFINES		= $(THREADS) \
		  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS	= $(INCDIRS:%=-I%) $(DEFINES)
CFLAGS		= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)
CXXFLAGS	=$(CXXDEBUG) $(CXXOPTIMISE) $(CXXANSI_CONFORMANCE)
 

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

soarchive:		$(SOARCHIVE)

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
			$(ARCHIVE) $(SOARCHIVE)

# Target to clean up everything that might not be wanted (should not need
# modifying).
clobber:		clean
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			-$(RM) *.o *.a .pure *.pure *.bak *~

# making the tag file
tags:
	$(ETAGS) -f $(ETAGSFILE) $(ETAGSFLAGS) $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)

allsources:		sources
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# Targets to create archives, dependencies and rules (modify as required).
$(ARCHIVE):		$(OBJECTS)
			$(AR) r $(ARCHIVE) $?

$(SOARCHIVE):		$(OBJECTS)
			$(CC) -G $(OBJECTS) -lC -o $(SOARCHIVE)

# Targets to create executables, just a dependencies (modify as required).
testWlzFacts:	testWlzFacts.o

# Target for Purified executables (modify as required).
purify:
testWlzFacts.pure:	testWlzFacts.o

$(OBJECTS):		%.o: %.c $(INCLUDES_ALL) 
ifeq ($(UNIXTYPE), SUNOS5)
			$(CC) -v  -c -o $@ $(CPPFLAGS) $(CFLAGS) $<
else
			$(CC) -c -o $@ $(CPPFLAGS) $(CFLAGS) $<
endif
