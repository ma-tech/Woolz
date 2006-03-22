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
# \brief        Makefile for the Woolz library.
# \ingroup      Wlz
# \todo         -
# \bug          None known.
##

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			= libWlz.a
SOARCHIVE		= $(ARCHIVE:%.a=%.so)
DLLARCHIVE		= $(ARCHIVE:%.a=cyg%.dll)
DLLIMPORTLIB		= $(ARCHIVE:%.a=%.dll.a)
PROTOTYPES		= WlzProto.h
DLLLIBDIRS		= $(DEPENDDIRS:%=-L%)
DLLLIBDEPS		= -lAlg -lAlc
# DLLLIBSEXT		= -L. -lg -lm -lc
DLLLIBSEXT		= 

# Names of executables to be built (modify as required).
EXECUTABLES		= 

# List of all 'C' source files (modify as required).
CSOURCES		= \
			Wlz2DContains.c \
			Wlz3DProjection.c \
			Wlz3DSection.c \
			Wlz3DSectionFromGeoModel.c \
			Wlz3DSectionSegmentObject.c \
			Wlz3DViewStructUtils.c \
			Wlz3DViewTransformObj.c \
			Wlz3DWarpMQ_S.c \
			WlzAffineTransform.c \
			WlzAffineTransformLSq.c \
			WlzArea.c \
			WlzArray.c \
                        WlzAssign.c \
			WlzBackground.c \
			WlzBasisFn.c \
			WlzBasisFnTransform.c \
			WlzBoundaryUtils.c \
			WlzBoundingBox.c \
			WlzBoundToObj.c \
			WlzCannyDeriche.c \
			WlzCbThreshold.c \
			WlzCCor.c \
			WlzCentreOfMass.c \
			WlzClipObjToBox.c \
			WlzCMeshTransform.c \
			WlzCompThresh.c \
			WlzConstruct3D.c \
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
			WlzDistMetric.c \
			WlzDistTransform.c \
			WlzDomainFill.c \
			WlzDomainUtils.c \
			WlzEmpty.c \
			WlzErosion.c \
			WlzExplode3D.c \
			WlzErosion4.c \
			WlzError.c \
			WlzFacts.c \
			WlzFillBlankPlanes.c \
			WlzFreeSpace.c \
			WlzGauss.c \
			WlzGaussNoise.c \
			WlzGeoModel.c \
			WlzGeoModelFilters.c \
			WlzGeoModelStats.c \
			WlzGeometry.c \
			WlzGeometryTrackUpAndDown_s.c \
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
			WlzGreySetRangeLut.c \
			WlzGreySetValue.c \
			WlzGreyStats.c \
			WlzGreyTemplate.c \
			WlzGreyTransfer.c \
			WlzGreyValue.c \
			WlzGreyValueMixing_s.c \
			WlzHasIntersect.c \
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
			WlzLBTDomain.c \
			WlzLabel.c \
			WlzLabel3d.c \
			WlzLaplacian.c \
			WlzLineArea.c \
			WlzMakeAffineTransform.c \
			WlzMakeCompound.c \
			WlzMakeIntervalValues.c \
			WlzMakeProperties.c \
			WlzMakeStructs.c \
			WlzMeshGen.c \
			WlzMeshTransform.c \
			WlzMeshUtils.c \
			WlzMatchICP.c \
			WlzMwrAngle.c \
			WlzNMSuppress.c \
			WlzObjToBoundary.c \
			WlzPoints.c \
			WlzPolarSample.c  \
			WlzPolyDecimate.c \
			WlzPolyEquispace.c \
			WlzPolyReverse.c \
			WlzPolySmooth.c \
			WlzPolyToObj.c \
			WlzPolyUtils.c \
			WlzPrinicipalAngle.c \
			WlzRank.c \
			WlzRaster.c \
			WlzReadObj.c \
			WlzRegCCor.c \
			WlzRegICP.c \
			WlzRGBAConvert.c \
			WlzRGBAGreyStats.c \
			WlzRGBAImageArithmetic.c \
			WlzRGBAModGradient.c \
			WlzRGBARange.c \
			WlzRGBAPixelUtils.c \
			WlzRGBAThreshold.c \
			WlzRsvFilter.c \
			WlzSampleObj.c \
			WlzSampleValuesAndCoords.c \
			WlzScalarArithmeticOp.c \
			WlzScalarBinaryOp.c \
			WlzScalarFeatures.c \
			WlzSepTrans.c \
			WlzSeqPar.c \
			WlzShadeCorrect.c \
			WlzShift.c \
			WlzSkeleton.c \
			WlzSnapFit.c \
			WlzSobel.c \
			WlzSplitObj.c \
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
			WlzVerticies.c \
			WlzVolume.c \
			WlzWindow.c \
			WlzWriteObj.c 
#			WlzXdr.c


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

# Source man documentation directory
SRC_MAN_DIR		= man

# Source html documentation directory
SRC_WWW_DIR		= html

# Destination html documentation directory
DST_WWW_DIR		= $(HGU_WWW_DIR)/src/Core/libWlz

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
DEFINES		= $(THREADS) \
		  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS	= $(INCDIRS:%=-I%) $(DEFINES)
CFLAGS		= $(CDEBUG) $(COPTIMISE) $(ANSI_CONFORMANCE)
CXXFLAGS	=$(CXXDEBUG) $(CXXOPTIMISE) $(CXXANSI_CONFORMANCE)
ifeq ($(UNIXTYPE), SUNOS5)
ANSI_CONFORMANCE	=
endif

# Default traget which builds everything (should not need modifying).
all:			archive executables includes documentation
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

# Target which installs everything after it has been built (should not need
# modifying).
install:		install_includes install_archive \
			install_executables install_documentation

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

# making the tag file
tags:			 $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)
			 $(CTAGS) $(CSOURCES) $(CCSOURCES) $(INCLUDES_ALL)

allsources:		sources
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif

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

# Targets to create executables, just a dependencies (modify as required).
testWlzFacts:	testWlzFacts.o

# Target for Purified executables (modify as required).
purify:
testWlzFacts.pure:	testWlzFacts.o

$(OBJECTS):		%.o: %.c $(INCLUDES_ALL) 
			$(CC) -c -o $@ $(CPPFLAGS) $(CFLAGS) $<



