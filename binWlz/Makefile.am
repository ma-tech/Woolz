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

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

LDADD 			= -L$(top_srcdir)/libWlz/.libs -lWlz \
			  -L$(top_srcdir)/libAlg/.libs -lAlg \
			  -L$(top_srcdir)/libAlc/.libs -lAlc \
			  -lm

bin_PROGRAMS		= \
			  Wlz3DWarpMQ \
			  WlzAffineTransformDiff \
			  WlzAffineTransformLSq \
			  WlzAffineTransformObj \
			  WlzAffineTransformVertices \
			  WlzArea \
			  WlzAutoCorrelate \
			  WlzBasisFnTransformObj \
			  WlzBasisFnTransformVertices \
			  WlzBoundaryToObj \
			  WlzBoundingBox \
			  WlzBuildContour \
			  WlzCannyDeriche \
			  WlzCentralityObj \
			  WlzCentreOfMass \
			  WlzChangeEMAPProperty \
			  WlzClipObjToBox \
			  WlzClosestVertices \
			  WlzCMeshDistance \
			  WlzCMeshGen \
			  WlzCMeshToSpatialDomain \
			  WlzCMeshTransformVtx \
			  WlzCMeshTransformObj \
			  WlzCompound \
			  WlzCompoundArrayToScalar \
			  WlzCompDispIncGrey \
			  WlzCompThresh \
			  WlzConstruct3D \
			  WlzContourFlipOrient \
			  WlzContourFromPoints \
			  WlzContourGeomFilter \
			  WlzContourObj \
			  WlzContourSpxStats \
			  WlzContourSzSelect \
			  WlzConvertPix \
			  WlzConvexPoly \
			  WlzCopyObj \
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
			  WlzFile \
			  WlzFilterNObjsValues \
			  WlzGauss \
			  WlzGaussNoise \
			  WlzGeometryTrackUpAndDown \
			  WlzGetSubSectionFromObj \
			  WlzGreyCrossing \
			  WlzGreyGradient \
			  WlzGreyInvert \
			  WlzGreyMask \
			  WlzGreyNormalise \
			  WlzGreyRange \
			  WlzGreySetIncValues \
			  WlzGreySetRange \
			  WlzGreySetValue \
			  WlzGreyStats \
			  WlzGreyTemplate \
			  WlzGreyTransfer \
			  WlzGreyValue \
			  WlzGreyValueMixing \
			  WlzGreyVarianceObj \
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
			  WlzIntersect \
			  WlzIntRescaleObj \
			  WlzIsEmpty \
			  WlzLabel \
			  WlzLaplacian \
			  WlzMakeEmpty \
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
			  WlzRGBCombine \
			  WlzRGBChannelRatio \
			  WlzRMSDiff \
			  WlzRsvFilterObj \
			  WlzSampleObj \
			  WlzScalarBinaryOp \
			  WlzScalarDivide \
			  WlzScalarFeatures \
			  WlzScalarFnObj \
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
			  WlzThreshold \
			  WlzTransformProduct \
			  WlzTransposeObj \
			  WlzTransposePlanes \
			  WlzUnion \
			  WlzVerifyObj \
			  WlzVolume \
			  WlzWindow

#			  WlzAutoTrackUpDown

Wlz3DWarpMQ_SOURCES			= Wlz3DWarpMQ.c
Wlz3DWarpMQ_LDADD			= $(LDADD)
Wlz3DWarpMQ_LDFLAGS			= $(AM_LFLAGS)

WlzAffineTransformDiff_SOURCES		= WlzAffineTransformDiff.c
WlzAffineTransformDiff_LDADD		= $(LDADD)
WlzAffineTransformDiff_LDFLAGS		= $(AM_LFLAGS)

WlzAffineTransformLSq_SOURCES		= WlzAffineTransformLSq.c
WlzAffineTransformLSq_LDADD		= $(LDADD)
WlzAffineTransformLSq_LDFLAGS		= $(AM_LFLAGS)

WlzAffineTransformVertices_SOURCES	= WlzAffineTransformVertices.c
WlzAffineTransformVertices_LDADD	= $(LDADD)
WlzAffineTransformVertices_LDFLAGS	= $(AM_LFLAGS)

WlzAffineTransformObj_SOURCES		= WlzAffineTransformObj.c
WlzAffineTransformObj_LDADD		= $(LDADD)
WlzAffineTransformObj_LDFLAGS		= $(AM_LFLAGS)

WlzArea_SOURCES				= WlzArea.c
WlzArea_LDADD				= $(LDADD)
WlzArea_LDFLAGS				= $(AM_LFLAGS)

WlzAutoCorrelate_SOURCES		= WlzAutoCorrelate.c
WlzAutoCorrelate_LDADD			= $(LDADD)
WlzAutoCorrelate_LDFLAGS		= $(AM_LFLAGS)

#WlzAutoTrackUpDown_SOURCES		= WlzAutoTrackUpDown.c
#WlzAutoTrackUpDown_LDADD		= $(LDADD)
#WlzAutoTrackUpDown_LDFLAGS		= $(AM_LFLAGS)

WlzBasisFnTransformObj_SOURCES		= WlzBasisFnTransformObj.c
WlzBasisFnTransformObj_LDADD		= $(LDADD)
WlzBasisFnTransformObj_LDFLAGS		= $(AM_LFLAGS)

WlzBasisFnTransformVertices_SOURCES	= WlzBasisFnTransformVertices.c
WlzBasisFnTransformVertices_LDADD	= $(LDADD)
WlzBasisFnTransformVertices_LDFLAGS	= $(AM_LFLAGS)

WlzBoundaryToObj_SOURCES		= WlzBoundaryToObj.c
WlzBoundaryToObj_LDADD			= $(LDADD)
WlzBoundaryToObj_LDFLAGS		= $(AM_LFLAGS)

WlzBoundingBox_SOURCES			= WlzBoundingBox.c
WlzBoundingBox_LDADD			= $(LDADD)
WlzBoundingBox_LDFLAGS			= $(AM_LFLAGS)

WlzBuildContour_SOURCES			= WlzBuildContour.c
WlzBuildContour_LDADD			= $(LDADD)
WlzBuildContour_LDFLAGS			= $(AM_LFLAGS)

WlzCannyDeriche_SOURCES			= WlzCannyDeriche.c
WlzCannyDeriche_LDADD			= $(LDADD)
WlzCannyDeriche_LDFLAGS			= $(AM_LFLAGS)

WlzCentralityObj_SOURCES		= WlzCentralityObj.c
WlzCentralityObj_LDADD			= $(LDADD)
WlzCentralityObj_LDFLAGS		= $(AM_LFLAGS)

WlzCentreOfMass_SOURCES			= WlzCentreOfMass.c
WlzCentreOfMass_LDADD			= $(LDADD)
WlzCentreOfMass_LDFLAGS			= $(AM_LFLAGS)

WlzChangeEMAPProperty_SOURCES		= WlzChangeEMAPProperty.c
WlzChangeEMAPProperty_LDADD		= $(LDADD)
WlzChangeEMAPProperty_LDFLAGS		= $(AM_LFLAGS)

WlzClipObjToBox_SOURCES			= WlzClipObjToBox.c
WlzClipObjToBox_LDADD			= $(LDADD)
WlzClipObjToBox_LDFLAGS			= $(AM_LFLAGS)

WlzClosestVertices_SOURCES		= WlzClosestVertices.c
WlzClosestVertices_LDADD		= $(LDADD)
WlzClosestVertices_LDFLAGS		= $(AM_LFLAGS)

WlzCMeshDistance_SOURCES		= WlzCMeshDistance.c
WlzCMeshDistance_LDADD			= $(LDADD)
WlzCMeshDistance_LDFLAGS		= $(AM_LFLAGS)

WlzCMeshGen_SOURCES			= WlzCMeshGen.c
WlzCMeshGen_LDADD			= $(LDADD)
WlzCMeshGen_LDFLAGS			= $(AM_LFLAGS)

WlzCMeshToSpatialDomain_SOURCES		= WlzCMeshToSpatialDomain.c
WlzCMeshToSpatialDomain_LDADD		= $(LDADD)
WlzCMeshToSpatialDomain_LDFLAGS		= $(AM_LFLAGS)

WlzCMeshTransformVtx_SOURCES		= WlzCMeshTransformVtx.c
WlzCMeshTransformVtx_LDADD		= $(LDADD)
WlzCMeshTransformVtx_LDFLAGS		= $(AM_LFLAGS)

WlzCMeshTransformObj_SOURCES		= WlzCMeshTransformObj.c
WlzCMeshTransformObj_LDADD		= $(LDADD)
WlzCMeshTransformObj_LDFLAGS		= $(AM_LFLAGS)

WlzCompDispIncGrey_SOURCES		= WlzCompDispIncGrey.c
WlzCompDispIncGrey_LDADD		= $(LDADD)
WlzCompDispIncGrey_LDFLAGS		= $(AM_LFLAGS)

WlzCompThresh_SOURCES			= WlzCompThresh.c
WlzCompThresh_LDADD			= $(LDADD)
WlzCompThresh_LDFLAGS			= $(AM_LFLAGS)

WlzCompound_SOURCES			= WlzCompound.c
WlzCompound_LDADD			= $(LDADD)
WlzCompound_LDFLAGS			= $(AM_LFLAGS)

WlzCompoundArrayToScalar_SOURCES	= WlzCompoundArrayToScalar.c
WlzCompoundArrayToScalar_LDADD		= $(LDADD)
WlzCompoundArrayToScalar_LDFLAGS	= $(AM_LFLAGS)

WlzConstruct3D_SOURCES			= WlzConstruct3D.c
WlzConstruct3D_LDADD			= $(LDADD)
WlzConstruct3D_LDFLAGS			= $(AM_LFLAGS)

WlzContourFlipOrient_SOURCES		= WlzContourFlipOrient.c
WlzContourFlipOrient_LDADD		= $(LDADD)
WlzContourFlipOrient_LDFLAGS		= $(AM_LFLAGS)

WlzContourFromPoints_SOURCES		= WlzContourFromPoints.c
WlzContourFromPoints_LDADD		= $(LDADD)
WlzContourFromPoints_LDFLAGS		= $(AM_LFLAGS)

WlzContourGeomFilter_SOURCES		= WlzContourGeomFilter.c
WlzContourGeomFilter_LDADD		= $(LDADD)
WlzContourGeomFilter_LDFLAGS		= $(AM_LFLAGS)

WlzContourObj_SOURCES			= WlzContourObj.c
WlzContourObj_LDADD			= $(LDADD)
WlzContourObj_LDFLAGS			= $(AM_LFLAGS)

WlzContourSzSelect_SOURCES		= WlzContourSzSelect.c
WlzContourSzSelect_LDADD		= $(LDADD)
WlzContourSzSelect_LDFLAGS		= $(AM_LFLAGS)

WlzContourSpxStats_SOURCES		= WlzContourSpxStats.c
WlzContourSpxStats_LDADD		= $(LDADD)
WlzContourSpxStats_LDFLAGS		= $(AM_LFLAGS)

WlzConvertPix_SOURCES			= WlzConvertPix.c
WlzConvertPix_LDADD			= $(LDADD)
WlzConvertPix_LDFLAGS			= $(AM_LFLAGS)

WlzConvexPoly_SOURCES			= WlzConvexPoly.c
WlzConvexPoly_LDADD			= $(LDADD)
WlzConvexPoly_LDFLAGS			= $(AM_LFLAGS)

WlzCopyObj_SOURCES			= WlzCopyObj.c
WlzCopyObj_LDADD			= $(LDADD)
WlzCopyObj_LDFLAGS			= $(AM_LFLAGS)

WlzCreateSpecialSE_SOURCES		= WlzCreateSpecialSE.c
WlzCreateSpecialSE_LDADD		= $(LDADD)
WlzCreateSpecialSE_LDFLAGS		= $(AM_LFLAGS)

WlzCrossCorValue_SOURCES		= WlzCrossCorValue.c
WlzCrossCorValue_LDADD			= $(LDADD)
WlzCrossCorValue_LDFLAGS		= $(AM_LFLAGS)

WlzCutObjToBox_SOURCES			= WlzCutObjToBox.c
WlzCutObjToBox_LDADD			= $(LDADD)
WlzCutObjToBox_LDFLAGS			= $(AM_LFLAGS)

WlzDiffDomain_SOURCES			= WlzDiffDomain.c
WlzDiffDomain_LDADD			= $(LDADD)
WlzDiffDomain_LDFLAGS			= $(AM_LFLAGS)

WlzDilation_SOURCES			= WlzDilation.c
WlzDilation_LDADD			= $(LDADD)
WlzDilation_LDFLAGS			= $(AM_LFLAGS)

WlzDistanceMetric_SOURCES		= WlzDistanceMetric.c
WlzDistanceMetric_LDADD			= $(LDADD)
WlzDistanceMetric_LDFLAGS		= $(AM_LFLAGS)

WlzDistanceTransform_SOURCES		= WlzDistanceTransform.c
WlzDistanceTransform_LDADD		= $(LDADD)
WlzDistanceTransform_LDFLAGS		= $(AM_LFLAGS)

WlzDomain_SOURCES			= WlzDomain.c
WlzDomain_LDADD				= $(LDADD)
WlzDomain_LDFLAGS			= $(AM_LFLAGS)

WlzDomainFill_SOURCES			= WlzDomainFill.c
WlzDomainFill_LDADD			= $(LDADD)
WlzDomainFill_LDFLAGS			= $(AM_LFLAGS)
				
WlzErosion_SOURCES			= WlzErosion.c
WlzErosion_LDADD			= $(LDADD)
WlzErosion_LDFLAGS			= $(AM_LFLAGS)

WlzExplode_SOURCES			= WlzExplode.c
WlzExplode_LDADD			= $(LDADD)
WlzExplode_LDFLAGS			= $(AM_LFLAGS)

WlzFacts_SOURCES			= WlzFacts.c
WlzFacts_LDADD				= $(LDADD)
WlzFacts_LDFLAGS			= $(AM_LFLAGS)

WlzFile_SOURCES				= WlzFile.c
WlzFile_LDADD				= $(LDADD)
WlzFile_LDFLAGS				= $(AM_LFLAGS)

WlzFilterNObjsValues_SOURCES		= WlzFilterNObjsValues.c
WlzFilterNObjsValues_LDADD		= $(LDADD)
WlzFilterNObjsValues_LDFLAGS		= $(AM_LFLAGS)

WlzGauss_SOURCES			= WlzGauss.c
WlzGauss_LDADD				= $(LDADD)
WlzGauss_LDFLAGS			= $(AM_LFLAGS)

WlzGaussNoise_SOURCES			= WlzGaussNoise.c
WlzGaussNoise_LDADD			= $(LDADD)
WlzGaussNoise_LDFLAGS			= $(AM_LFLAGS)

WlzGeometryTrackUpAndDown_SOURCES	= WlzGeometryTrackUpAndDown.c
WlzGeometryTrackUpAndDown_LDADD		= $(LDADD)
WlzGeometryTrackUpAndDown_LDFLAGS	= $(AM_LFLAGS)

WlzGetSubSectionFromObj_SOURCES		= WlzGetSubSectionFromObj.c
WlzGetSubSectionFromObj_LDADD		= $(LDADD)
WlzGetSubSectionFromObj_LDFLAGS		= $(AM_LFLAGS)

WlzGreyCrossing_SOURCES			= WlzGreyCrossing.c
WlzGreyCrossing_LDADD			= $(LDADD)
WlzGreyCrossing_LDFLAGS			= $(AM_LFLAGS)

WlzGreyGradient_SOURCES			= WlzGreyGradient.c
WlzGreyGradient_LDADD			= $(LDADD)
WlzGreyGradient_LDFLAGS			= $(AM_LFLAGS)

WlzGreyInvert_SOURCES			= WlzGreyInvert.c
WlzGreyInvert_LDADD			= $(LDADD)
WlzGreyInvert_LDFLAGS			= $(AM_LFLAGS)

WlzGreyMask_SOURCES			= WlzGreyMask.c
WlzGreyMask_LDADD			= $(LDADD)
WlzGreyMask_LDFLAGS			= $(AM_LFLAGS)

WlzGreyNormalise_SOURCES		= WlzGreyNormalise.c
WlzGreyNormalise_LDADD			= $(LDADD)
WlzGreyNormalise_LDFLAGS		= $(AM_LFLAGS)

WlzGreyRange_SOURCES			= WlzGreyRange.c
WlzGreyRange_LDADD			= $(LDADD)
WlzGreyRange_LDFLAGS			= $(AM_LFLAGS)

WlzGreySetIncValues_SOURCES		= WlzGreySetIncValues.c
WlzGreySetIncValues_LDADD		= $(LDADD)
WlzGreySetIncValues_LDFLAGS		= $(AM_LFLAGS)

WlzGreySetRange_SOURCES			= WlzGreySetRange.c
WlzGreySetRange_LDADD			= $(LDADD)
WlzGreySetRange_LDFLAGS			= $(AM_LFLAGS)

WlzGreySetValue_SOURCES			= WlzGreySetValue.c
WlzGreySetValue_LDADD			= $(LDADD)
WlzGreySetValue_LDFLAGS			= $(AM_LFLAGS)

WlzGreyStats_SOURCES			= WlzGreyStats.c
WlzGreyStats_LDADD			= $(LDADD)
WlzGreyStats_LDFLAGS			= $(AM_LFLAGS)

WlzGreyTemplate_SOURCES			= WlzGreyTemplate.c
WlzGreyTemplate_LDADD			= $(LDADD)
WlzGreyTemplate_LDFLAGS			= $(AM_LFLAGS)

WlzGreyTransfer_SOURCES			= WlzGreyTransfer.c
WlzGreyTransfer_LDADD			= $(LDADD)
WlzGreyTransfer_LDFLAGS			= $(AM_LFLAGS)

WlzGreyValue_SOURCES			= WlzGreyValue.c
WlzGreyValue_LDADD			= $(LDADD)
WlzGreyValue_LDFLAGS			= $(AM_LFLAGS)

WlzGreyValueMixing_SOURCES		= WlzGreyValueMixing.c
WlzGreyValueMixing_LDADD		= $(LDADD)
WlzGreyValueMixing_LDFLAGS		= $(AM_LFLAGS)

WlzGreyVarianceObj_SOURCES		= WlzGreyVarianceObj.c
WlzGreyVarianceObj_LDADD		= $(LDADD)
WlzGreyVarianceObj_LDFLAGS		= $(AM_LFLAGS)

WlzHasIntersect_SOURCES			= WlzHasIntersect.c
WlzHasIntersect_LDADD			= $(LDADD)
WlzHasIntersect_LDFLAGS			= $(AM_LFLAGS)

WlzHistogramData_SOURCES		= WlzHistogramData.c
WlzHistogramData_LDADD			= $(LDADD)
WlzHistogramData_LDFLAGS		= $(AM_LFLAGS)

WlzHistogramEqualiseObj_SOURCES		= WlzHistogramEqualiseObj.c
WlzHistogramEqualiseObj_LDADD		= $(LDADD)
WlzHistogramEqualiseObj_LDFLAGS		= $(AM_LFLAGS)

WlzHistogramFindPeaks_SOURCES		= WlzHistogramFindPeaks.c
WlzHistogramFindPeaks_LDADD		= $(LDADD)
WlzHistogramFindPeaks_LDFLAGS		= $(AM_LFLAGS)

WlzHistogramFitPeaks_SOURCES		= WlzHistogramFitPeaks.c
WlzHistogramFitPeaks_LDADD		= $(LDADD)
WlzHistogramFitPeaks_LDFLAGS		= $(AM_LFLAGS)

WlzHistogramGauss_SOURCES		= WlzHistogramGauss.c
WlzHistogramGauss_LDADD			= $(LDADD)
WlzHistogramGauss_LDFLAGS		= $(AM_LFLAGS)

WlzHistogramMatchObj_SOURCES		= WlzHistogramMatchObj.c
WlzHistogramMatchObj_LDADD		= $(LDADD)
WlzHistogramMatchObj_LDFLAGS		= $(AM_LFLAGS)

WlzHistogramObj_SOURCES			= WlzHistogramObj.c
WlzHistogramObj_LDADD			= $(LDADD)
WlzHistogramObj_LDFLAGS			= $(AM_LFLAGS)

WlzHistogramRebin_SOURCES		= WlzHistogramRebin.c
WlzHistogramRebin_LDADD			= $(LDADD)
WlzHistogramRebin_LDFLAGS		= $(AM_LFLAGS)

WlzHistogramSmooth_SOURCES		= WlzHistogramSmooth.c
WlzHistogramSmooth_LDADD		= $(LDADD)
WlzHistogramSmooth_LDFLAGS		= $(AM_LFLAGS)

WlzImageArithmetic_SOURCES		= WlzImageArithmetic.c
WlzImageArithmetic_LDADD		= $(LDADD)
WlzImageArithmetic_LDFLAGS		= $(AM_LFLAGS)

WlzInsideDomain_SOURCES			= WlzInsideDomain.c
WlzInsideDomain_LDADD			= $(LDADD)
WlzInsideDomain_LDFLAGS			= $(AM_LFLAGS)

WlzIntRescaleObj_SOURCES		= WlzIntRescaleObj.c
WlzIntRescaleObj_LDADD			= $(LDADD)
WlzIntRescaleObj_LDFLAGS		= $(AM_LFLAGS)

WlzIsEmpty_SOURCES			= WlzIsEmpty.c
WlzIsEmpty_LDADD			= $(LDADD)
WlzIsEmpty_LDFLAGS			= $(AM_LFLAGS)

WlzIntersect_SOURCES			= WlzIntersect.c
WlzIntersect_LDADD			= $(LDADD)
WlzIntersect_LDFLAGS			= $(AM_LFLAGS)

WlzLabel_SOURCES			= WlzLabel.c
WlzLabel_LDADD				= $(LDADD)
WlzLabel_LDFLAGS			= $(AM_LFLAGS)

WlzLaplacian_SOURCES			= WlzLaplacian.c
WlzLaplacian_LDADD			= $(LDADD)
WlzLaplacian_LDFLAGS			= $(AM_LFLAGS)

WlzMakeEmpty_SOURCES			= WlzMakeEmpty.c
WlzMakeEmpty_LDADD			= $(LDADD)
WlzMakeEmpty_LDFLAGS			= $(AM_LFLAGS)

WlzMarkersToDomain_SOURCES		= WlzMarkersToDomain.c
WlzMarkersToDomain_LDADD		= $(LDADD)
WlzMarkersToDomain_LDFLAGS		= $(AM_LFLAGS)

WlzMatchICPObj_SOURCES			= WlzMatchICPObj.c
WlzMatchICPObj_LDADD			= $(LDADD)
WlzMatchICPObj_LDFLAGS			= $(AM_LFLAGS)

WlzMeshTransformObj_SOURCES		= WlzMeshTransformObj.c
WlzMeshTransformObj_LDADD		= $(LDADD)
WlzMeshTransformObj_LDFLAGS		= $(AM_LFLAGS)

WlzObjToBoundary_SOURCES		= WlzObjToBoundary.c
WlzObjToBoundary_LDADD			= $(LDADD)
WlzObjToBoundary_LDFLAGS		= $(AM_LFLAGS)

WlzPolarSample_SOURCES			= WlzPolarSample.c
WlzPolarSample_LDADD			= $(LDADD)
WlzPolarSample_LDFLAGS			= $(AM_LFLAGS)

WlzPolygonToObj_SOURCES			= WlzPolygonToObj.c
WlzPolygonToObj_LDADD			= $(LDADD)
WlzPolygonToObj_LDFLAGS			= $(AM_LFLAGS)

WlzPrinicipalAngle_SOURCES		= WlzPrinicipalAngle.c
WlzPrinicipalAngle_LDADD		= $(LDADD)
WlzPrinicipalAngle_LDFLAGS		= $(AM_LFLAGS)

WlzRandomAffineTransform_SOURCES	= WlzRandomAffineTransform.c
WlzRandomAffineTransform_LDADD		= $(LDADD)
WlzRandomAffineTransform_LDFLAGS	= $(AM_LFLAGS)

WlzRankObj_SOURCES			= WlzRankObj.c
WlzRankObj_LDADD			= $(LDADD)
WlzRankObj_LDFLAGS			= $(AM_LFLAGS)

WlzRasterObj_SOURCES			= WlzRasterObj.c
WlzRasterObj_LDADD			= $(LDADD)
WlzRasterObj_LDFLAGS			= $(AM_LFLAGS)

WlzRBFContourObj_SOURCES		= WlzRBFContourObj.c
WlzRBFContourObj_LDADD			= $(LDADD)
WlzRBFContourObj_LDFLAGS		= $(AM_LFLAGS)

WlzRegisterCCor_SOURCES			= WlzRegisterCCor.c
WlzRegisterCCor_LDADD			= $(LDADD)
WlzRegisterCCor_LDFLAGS			= $(AM_LFLAGS)

WlzRegisterICP_SOURCES			= WlzRegisterICP.c
WlzRegisterICP_LDADD			= $(LDADD)
WlzRegisterICP_LDFLAGS			= $(AM_LFLAGS)

WlzRegisterICPWSD_SOURCES		= WlzRegisterICPWSD.c
WlzRegisterICPWSD_LDADD			= $(LDADD)
WlzRegisterICPWSD_LDFLAGS		= $(AM_LFLAGS)

WlzRGBAConvert_SOURCES			= WlzRGBAConvert.c
WlzRGBAConvert_LDADD			= $(LDADD)
WlzRGBAConvert_LDFLAGS			= $(AM_LFLAGS)

WlzRGBChannelRatio_SOURCES		= WlzRGBChannelRatio.c
WlzRGBChannelRatio_LDADD		= $(LDADD)
WlzRGBChannelRatio_LDFLAGS		= $(AM_LFLAGS)

WlzRGBCombine_SOURCES			= WlzRGBCombine.c
WlzRGBCombine_LDADD			= $(LDADD)
WlzRGBCombine_LDFLAGS			= $(AM_LFLAGS)

WlzRMSDiff_SOURCES			= WlzRMSDiff.c
WlzRMSDiff_LDADD			= $(LDADD)
WlzRMSDiff_LDFLAGS			= $(AM_LFLAGS)

WlzRsvFilterObj_SOURCES			= WlzRsvFilterObj.c
WlzRsvFilterObj_LDADD			= $(LDADD)
WlzRsvFilterObj_LDFLAGS			= $(AM_LFLAGS)

WlzSampleObj_SOURCES			= WlzSampleObj.c
WlzSampleObj_LDADD			= $(LDADD)
WlzSampleObj_LDFLAGS			= $(AM_LFLAGS)

WlzScalarBinaryOp_SOURCES		= WlzScalarBinaryOp.c
WlzScalarBinaryOp_LDADD			= $(LDADD)
WlzScalarBinaryOp_LDFLAGS		= $(AM_LFLAGS)

WlzScalarDivide_SOURCES			= WlzScalarDivide.c
WlzScalarDivide_LDADD			= $(LDADD)
WlzScalarDivide_LDFLAGS			= $(AM_LFLAGS)

WlzScalarFeatures_SOURCES		= WlzScalarFeatures.c
WlzScalarFeatures_LDADD			= $(LDADD)
WlzScalarFeatures_LDFLAGS		= $(AM_LFLAGS)

WlzScalarFnObj_SOURCES			= WlzScalarFnObj.c
WlzScalarFnObj_LDADD			= $(LDADD)
WlzScalarFnObj_LDFLAGS			= $(AM_LFLAGS)

WlzSelect1InN_SOURCES			= WlzSelect1InN.c
WlzSelect1InN_LDADD			= $(LDADD)
WlzSelect1InN_LDFLAGS			= $(AM_LFLAGS)

WlzSetBackground_SOURCES		= WlzSetBackground.c
WlzSetBackground_LDADD			= $(LDADD)
WlzSetBackground_LDFLAGS		= $(AM_LFLAGS)

WlzSetVoxelSize_SOURCES			= WlzSetVoxelSize.c
WlzSetVoxelSize_LDADD			= $(LDADD)
WlzSetVoxelSize_LDFLAGS			= $(AM_LFLAGS)

WlzShadeCorrect_SOURCES			= WlzShadeCorrect.c
WlzShadeCorrect_LDADD			= $(LDADD)
WlzShadeCorrect_LDFLAGS			= $(AM_LFLAGS)

WlzShiftObj_SOURCES			= WlzShiftObj.c
WlzShiftObj_LDADD			= $(LDADD)
WlzShiftObj_LDFLAGS			= $(AM_LFLAGS)

WlzSkeleton_SOURCES			= WlzSkeleton.c
WlzSkeleton_LDADD			= $(LDADD)
WlzSkeleton_LDFLAGS			= $(AM_LFLAGS)

WlzSnapFitObjs_SOURCES			= WlzSnapFitObjs.c
WlzSnapFitObjs_LDADD			= $(LDADD)
WlzSnapFitObjs_LDFLAGS			= $(AM_LFLAGS)

WlzSobel_SOURCES			= WlzSobel.c
WlzSobel_LDADD				= $(LDADD)
WlzSobel_LDFLAGS			= $(AM_LFLAGS)

WlzStripProperties_SOURCES		= WlzStripProperties.c
WlzStripProperties_LDADD		= $(LDADD)
WlzStripProperties_LDFLAGS		= $(AM_LFLAGS)

WlzStructDilation_SOURCES		= WlzStructDilation.c
WlzStructDilation_LDADD			= $(LDADD)
WlzStructDilation_LDFLAGS		= $(AM_LFLAGS)

WlzStructErosion_SOURCES		= WlzStructErosion.c
WlzStructErosion_LDADD			= $(LDADD)
WlzStructErosion_LDFLAGS		= $(AM_LFLAGS)

WlzTransformProduct_SOURCES		= WlzTransformProduct.c
WlzTransformProduct_LDADD		= $(LDADD)
WlzTransformProduct_LDFLAGS		= $(AM_LFLAGS)

WlzTransposeObj_SOURCES			= WlzTransposeObj.c
WlzTransposeObj_LDADD			= $(LDADD)
WlzTransposeObj_LDFLAGS			= $(AM_LFLAGS)

WlzTransposePlanes_SOURCES		= WlzTransposePlanes.c
WlzTransposePlanes_LDADD		= $(LDADD)
WlzTransposePlanes_LDFLAGS		= $(AM_LFLAGS)

WlzThreshold_SOURCES			= WlzThreshold.c
WlzThreshold_LDADD			= $(LDADD)
WlzThreshold_LDFLAGS			= $(AM_LFLAGS)

WlzUnion_SOURCES			= WlzUnion.c
WlzUnion_LDADD				= $(LDADD)
WlzUnion_LDFLAGS			= $(AM_LFLAGS)

WlzVerifyObj_SOURCES			= WlzVerifyObj.c
WlzVerifyObj_LDADD			= $(LDADD)
WlzVerifyObj_LDFLAGS			= $(AM_LFLAGS)

WlzVolume_SOURCES			= WlzVolume.c
WlzVolume_LDADD				= $(LDADD)
WlzVolume_LDFLAGS			= $(AM_LFLAGS)

WlzWindow_SOURCES			= WlzWindow.c
WlzWindow_LDADD				= $(LDADD)
WlzWindow_LDFLAGS			= $(AM_LFLAGS)

