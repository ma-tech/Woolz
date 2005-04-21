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
# 13-01-04 Rao  add WlzBasisFnTransformVertices.c
# 23-04-03 Rao  Add WlzGreyValueMixing.c#
# 28-01-03 Rao  Add WlzGeometryTrackUpAndDown.c
# 18-09-02 Jianguo add Wlz3DWarpMQ
# 14-03-01 bill Add WlzRasterObj.
# 08-03-01 bill Add WlzContourSzSelect.
# 06-12-00 bill Add WlzRegisterICP.
# 16-08-00 bill Add WlzContourObj.
# 07-02-00 bill	Add WlzHistogramFitPeaks.
# 02-02-00 bill	Add WlzHistogramFindPeaks and WlzHistogramGauss.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			=

# Names of executables to be built (modify as required).
EXECUTABLES		= \
			WlzAffineTransformDiff \
			WlzAffineTransformLSq \
			WlzAffineTransformVertices \
			WlzAffineTransformObj \
			WlzArea \
			WlzBasisFnTransformObj \
			WlzBasisFnTransformVertices \
			WlzBoundaryToObj \
			WlzBoundingBox \
			WlzBuildContour \
			WlzCannyDeriche \
			WlzCentreOfMass \
			WlzChangeEMAPProperty \
			WlzClipObjToBox \
			WlzClosestVertices \
			WlzCompThresh \
			WlzCompound \
			WlzConstruct3D \
			WlzContourFlipOrient \
			WlzContourFromPoints \
			WlzContourGeomFilter \
			WlzContourObj \
			WlzContourSzSelect \
			WlzContourSpxStats \
			WlzConvertPix \
			WlzConvexPoly \
			WlzCreateSpecialSE \
			WlzCrossCorValue \
			WlzCutObjToBox \
			WlzDiffDomain \
			WlzDilation \
			WlzDistanceMetric \
			WlzDistanceTransform \
			WlzDomain \
			WlzDomainFill \
			WlzErosion \
			WlzExplode \
			WlzFacts \
			WlzGauss \
			WlzGaussNoise \
			WlzGeometryTrackUpAndDown \
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
			WlzGreyTransfer \
			WlzGreyValue \
	                WlzGreyValueMixing \
		        WlzHasIntersect \
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
			WlzInsideDomain \
			WlzIntRescaleObj \
			WlzIsEmpty \
			WlzIntersect \
			WlzLabel \
			WlzLaplacian \
			WlzMarkersToDomain \
			WlzMatchICPObj \
			WlzMeshTransformObj \
			WlzObjToBoundary \
			WlzPolarSample \
			WlzPolygonToObj \
			WlzPrinicipalAngle \
			WlzRandomAffineTransform \
			WlzRankObj \
			WlzRasterObj \
			WlzRBFContourObj \
			WlzRegisterCCor \
			WlzRegisterICP \
			WlzRegisterICPWSD \
			WlzRGBAConvert \
			WlzRMSDiff \
			WlzRsvFilterObj \
			WlzSampleObj \
			WlzScalarFeatures \
			WlzSelect1InN \
			WlzSetBackground \
			WlzSetVoxelSize \
			WlzShadeCorrect \
			WlzShiftObj \
			WlzSkeleton \
			WlzSnapFitObjs \
			WlzSobel \
			WlzStripProperties \
			WlzStructDilation \
			WlzStructErosion \
			WlzTransformProduct \
			WlzTransposeObj \
			WlzTransposePlanes \
			WlzThreshold \
			WlzUnion \
			WlzVerifyObj \
			WlzVolume \
			WlzWindow \
			Wlz3DWarpMQ

#			WlzAutoTrackUpDown



# List of all 'C' source files (modify as required).
CSOURCES		= \
			Wlz3DWarpMQ.c \
			MADomainTest.c \
			WlzAffineTransformDiff.c \
			WlzAffineTransformLSq.c \
			WlzAffineTransformVertices.c \
			WlzAffineTransformObj.c \
			WlzArea.c \
			WlzBasisFnTransformObj.c \
			WlzBasisFnTransformVertices.c \
			WlzBoundaryToObj.c \
			WlzBoundingBox.c \
			WlzBuildContour.c \
			WlzCannyDeriche.c \
			WlzCentreOfMass.c \
			WlzChangeEMAPProperty.c \
			WlzClipObjToBox.c \
			WlzClosestVertices.c \
			WlzCompound.c \
			WlzCompThresh.c \
			WlzConstruct3D.c \
			WlzContourFlipOrient.c \
			WlzContourFromPoints.c \
			WlzContourGeomFilter.c \
			WlzContourObj.c \
			WlzContourSpxStats.c \
			WlzContourSzSelect.c \
			WlzConvertPix.c \
			WlzConvexPoly.c \
			WlzCreateSpecialSE.c \
			WlzCrossCorValue.c \
			WlzCutObjToBox.c \
			WlzDiffDomain.c \
			WlzDilation.c \
			WlzDistanceMetric.c \
			WlzDistanceTransform.c \
			WlzDomain.c \
			WlzDomainFill.c \
			AlgDPSearch.c \
			WlzErosion.c \
			WlzExplode.c \
			WlzExplode3D.c \
			WlzFacts.c \
			WlzGauss.c \
			WlzGaussNoise.c \
			WlzGeometryTrackUpAndDown.c \
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
			WlzGreyTransfer.c \
			WlzGreyValue.c \
			WlzGreyValueMixing.c \
			WlzHasIntersect.c \
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
			WlzInsideDomain.c \
			WlzIntRescaleObj.c \
			WlzIsEmpty.c \
			WlzIntersect.c \
			WlzLabel.c \
			WlzLaplacian.c \
			WlzMarkersToDomain.c \
			WlzMatchICPObj.c \
			WlzMeshTransformObj.c \
			WlzObjToBoundary.c \
			WlzPolarSample.c \
			WlzPolygonToObj.c \
			WlzPrinicipalAngle.c \
			WlzRandomAffineTransform.c \
			WlzRankObj.c \
			WlzRasterObj.c \
			WlzRBFContourObj.c \
			WlzRegisterCCor.c \
			WlzRegisterICP.c \
			WlzRegisterICPWSD.c \
			WlzRGBAConvert.c \
			WlzRMSDiff.c \
			WlzRsvFilterObj.c \
			WlzSampleObj.c \
			WlzScalarFeatures.c \
			WlzSelect1InN.c \
			WlzSetBackground.c \
			WlzSetVoxelSize.c \
			WlzShadeCorrect.c \
			WlzShiftObj.c \
			WlzSkeleton.c \
			WlzSnapFitObjs.c \
			WlzSobel.c \
			WlzStripProperties.c \
			WlzStructDilation.c \
			WlzStructErosion.c \
			WlzTransformProduct.c \
			WlzTransposeObj.c \
			WlzTransposePlanes.c \
			WlzThreshold.c \
			WlzUnion.c \
			WlzVerifyObj.c \
			WlzVolume.c \
			WlzWindow.c

#			WlzAutoTrackUpDown.c

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
			  WlzTransposePlanes.1 \
			  WlzWindow.1

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  . \
			  ../libAlc \
			  ../libAlg \
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

LIBRARIES		= Wlz Alg Alc \
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
#DEFINES			= -DUSE_OSTORE
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
install_executables:	$(EXECUTABLES)
ifneq ($(strip $(SUBSYSTEMS)),)
			$(subsystems)
endif
			$(install_executables)

# Target which installs release executables (should not need modifying).
install_release:	$(EXECUTABLES)
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
			$(CC) -c -o $@ $(CPPFLAGS) $(CFLAGS) $<
