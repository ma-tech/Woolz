##
# \file         Makefile.am
# \author       Bill Hill
# \date         November 2007
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
#
# \brief        Process this file with automake to produce Makefile.in
# \todo         -
# \bug          None known.
##

AM_CPPFLAGS		= \
			  -I. \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

lib_LTLIBRARIES 	= libWlz.la

libWlz_la_LDFLAGS 	= -release @PACKAGE_VERSION@

libWlz_la_SOURCES 	= \
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
			  WlzAutoCor.c \
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
			  WlzCMeshUtils.c \
			  WlzCompThresh.c \
			  WlzConstruct3D.c \
			  WlzContour.c \
			  WlzConvertPix.c \
			  WlzConvexHull.c \
			  WlzConvolve.c \
			  WlzCopy.c \
			  WlzCutObjToBox.c \
			  WlzDebug.c \
			  WlzDiffDomain3d.c \
			  WlzDiffDomain.c \
			  WlzDilation.c \
			  WlzDistMetric.c \
			  WlzDistTransform.c \
			  WlzDomainFill.c \
			  WlzDomainUtils.c \
			  WlzEmpty.c \
			  WlzErosion4.c \
			  WlzErosion.c \
			  WlzError.c \
			  WlzExplode3D.c \
			  WlzFacts.c \
			  WlzFillBlankPlanes.c \
			  WlzFilterNObjValues.c \
			  WlzFreeSpace.c \
			  WlzGauss.c \
			  WlzGaussNoise.c \
			  WlzGeometry.c \
			  WlzGeometryTrackUpAndDown_s.c \
			  WlzGeoModel.c \
			  WlzGeoModelFilters.c \
			  WlzGeoModelStats.c \
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
			  WlzGreyVariance.c \
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
			  WlzLabel3d.c \
			  WlzLabel.c \
			  WlzLaplacian.c \
			  WlzLBTDomain.c \
			  WlzLineArea.c \
			  WlzMakeAffineTransform.c \
			  WlzMakeCompound.c \
			  WlzMakeIntervalValues.c \
			  WlzMakeProperties.c \
			  WlzMakeStructs.c \
			  WlzMatchICP.c \
			  WlzMeshGen.c \
			  WlzMeshTransform.c \
			  WlzMeshUtils.c \
			  WlzMwrAngle.c \
			  WlzNMSuppress.c \
			  WlzObjToBoundary.c \
			  WlzPoints.c \
			  WlzPolarSample.c \
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
			  WlzRGBAScalarBinaryOp.c \
			  WlzRGBAPixelUtils.c \
			  WlzRGBARange.c \
			  WlzRGBAThreshold.c \
			  WlzRGBChanRatio.c \
			  WlzRsvFilter.c \
			  WlzSampleObj.c \
			  WlzSampleValuesAndCoords.c \
			  WlzScalarArithmeticOp.c \
			  WlzScalarBinaryOp.c \
			  WlzScalarFeatures.c \
			  WlzScalarFn.c \
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
			  WlzTransform.c \
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

include_HEADERS 	= \
			  Wlz.h \
			  WlzError.h \
			  WlzDebug.h \
			  WlzProto.h \
			  WlzMacro.h \
			  WlzType.h
