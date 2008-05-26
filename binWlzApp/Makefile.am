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
			  -I$(top_srcdir)/libReconstruct \
			  -I$(top_srcdir)/libWlzExtFF \
			  -I$(top_srcdir)/libhguDlpList \
			  -I$(top_srcdir)/libbibfile \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
			  $(all_includes)

LDADD 			= \
			  -L$(top_srcdir)/libReconstruct/.libs -lReconstruct \
			  -L$(top_srcdir)/libWlzExtFF/.libs -lWlzExtFF \
			  -L$(top_srcdir)/libhguDlpList/.libs -lhguDlpList \
			  -L$(top_srcdir)/libbibfile/.libs -lbibfile \
			  -L$(top_srcdir)/libWlz/.libs -lWlz \
			  -L$(top_srcdir)/libAlg/.libs -lAlg \
			  -L$(top_srcdir)/libAlc/.libs -lAlc \
			  -ljpeg \
			  -ltiff \
			  -lm

bin_PROGRAMS		=  \
			  Wlz3DGetProjection \
			  Wlz3DGetSection \
			  Wlz3DViewTransformBitmap \
			  Wlz3DViewTransformObj \
			  WlzBibToMeshTrans \
			  WlzBibTransAbsRel \
			  WlzCelCmpFeatures \
			  WlzClassMap \
			  WlzCreateRandomDomain \
			  WlzDisplacement \
			  WlzDisplacementToVtk \
			  WlzDomainMatch \
			  WlzDomainMatchMatrix \
			  WlzEMAPDomainTransform \
			  WlzEuclidDistanceMap \
			  WlzExtractTransform \
			  WlzFixDomains \
			  WlzFixedPlaneAlign \
			  WlzIndexBoundaries \
			  WlzInterpolateSectionViewerBibFiles \
			  WlzMatchICPPlane \
			  WlzObjCompareSpecial_01 \
			  WlzOccupancy \
			  WlzPatchObjRegister \
			  WlzPosterProcessForStackupKeep3StraitLines \
			  WlzProductTransform \
			  WlzRawToWlz \
			  WlzReconstruct \
			  WlzSizeSelect \
			  WlzSplitImage \
			  WlzSplitMontage \
			  WlzTensor \
			  WlzTensorToVtk \
			  WlzTiffStackToPatches \
			  WlzTiffStackToShade \
			  WlzVtxDistance \
			  WlzVtxToObj

Wlz3DGetProjection_SOURCES	= Wlz3DGetProjection.c
Wlz3DGetProjection_LDADD	= $(LDADD)
Wlz3DGetProjection_LDFLAGS	= $(AM_LFLAGS)

Wlz3DGetSection_SOURCES		= Wlz3DGetSection.c
Wlz3DGetSection_LDADD		= $(LDADD)
Wlz3DGetSection_LDFLAGS		= $(AM_LFLAGS)

Wlz3DViewTransformBitmap_SOURCES = Wlz3DViewTransformBitmap.c
Wlz3DViewTransformBitmap_LDADD	= $(LDADD)
Wlz3DViewTransformBitmap_LDFLAGS = $(AM_LFLAGS)

Wlz3DViewTransformObj_SOURCES	= Wlz3DViewTransformObj.c
Wlz3DViewTransformObj_LDADD	= $(LDADD)
Wlz3DViewTransformObj_LDFLAGS	= $(AM_LFLAGS)

WlzBibToMeshTrans_SOURCES	= WlzBibToMeshTrans.c
WlzBibToMeshTrans_LDADD		= $(LDADD)
WlzBibToMeshTrans_LDFLAGS	= $(AM_LFLAGS)

WlzBibTransAbsRel_SOURCES	= WlzBibTransAbsRel.c
WlzBibTransAbsRel_LDADD		= $(LDADD)
WlzBibTransAbsRel_LDFLAGS	= $(AM_LFLAGS)

WlzCelCmpFeatures_SOURCES	= WlzCelCmpFeatures.c
WlzCelCmpFeatures_LDADD		= $(LDADD)
WlzCelCmpFeatures_LDFLAGS	= $(AM_LFLAGS)

WlzClassMap_SOURCES		= WlzClassMap.c
WlzClassMap_LDADD		= $(LDADD)
WlzClassMap_LDFLAGS		= $(AM_LFLAGS)

WlzCreateRandomDomain_SOURCES	= WlzCreateRandomDomain.c
WlzCreateRandomDomain_LDADD	= $(LDADD)
WlzCreateRandomDomain_LDFLAGS	= $(AM_LFLAGS)

WlzDisplacement_SOURCES		= WlzDisplacement.c
WlzDisplacement_LDADD		= $(LDADD)
WlzDisplacement_LDFLAGS		= $(AM_LFLAGS)

WlzDisplacementToVtk_SOURCES	= WlzDisplacementToVtk.c
WlzDisplacementToVtk_LDADD	= $(LDADD)
WlzDisplacementToVtk_LDFLAGS	= $(AM_LFLAGS)

WlzDomainMatch_SOURCES		= WlzDomainMatch.c
WlzDomainMatch_LDADD		= $(LDADD)
WlzDomainMatch_LDFLAGS		= $(AM_LFLAGS)

WlzDomainMatchMatrix_SOURCES	= WlzDomainMatchMatrix.c
WlzDomainMatchMatrix_LDADD	= $(LDADD)
WlzDomainMatchMatrix_LDFLAGS	= $(AM_LFLAGS)

WlzEMAPDomainTransform_SOURCES	= WlzEMAPDomainTransform.c \
				  WlzEMAPTransformUtils.c \
				  WlzEMAP.h
WlzEMAPDomainTransform_LDADD	= $(LDADD)
WlzEMAPDomainTransform_LDFLAGS	= $(AM_LFLAGS)

WlzEuclidDistanceMap_SOURCES	= WlzEuclidDistanceMap.c
WlzEuclidDistanceMap_LDADD	= $(LDADD)
WlzEuclidDistanceMap_LDFLAGS	= $(AM_LFLAGS)

WlzExtractTransform_SOURCES	= WlzExtractTransform.c
WlzExtractTransform_LDADD	= $(LDADD)
WlzExtractTransform_LDFLAGS	= $(AM_LFLAGS)

WlzFixDomains_SOURCES		= WlzFixDomains.c
WlzFixDomains_LDADD		= $(LDADD)
WlzFixDomains_LDFLAGS		= $(AM_LFLAGS)

WlzFixedPlaneAlign_SOURCES	= WlzFixedPlaneAlign.c
WlzFixedPlaneAlign_LDADD	= $(LDADD)
WlzFixedPlaneAlign_LDFLAGS	= $(AM_LFLAGS)

WlzIndexBoundaries_SOURCES	= WlzIndexBoundaries.c
WlzIndexBoundaries_LDADD	= $(LDADD)
WlzIndexBoundaries_LDFLAGS	= $(AM_LFLAGS)

WlzInterpolateSectionViewerBibFiles_SOURCES = \
                                          WlzInterpolateSectionViewerBibFiles.c
WlzInterpolateSectionViewerBibFiles_LDADD = $(LDADD)
WlzInterpolateSectionViewerBibFiles_LDFLAGS = $(AM_LFLAGS)

WlzMatchICPPlane_SOURCES	= WlzMatchICPPlane.c
WlzMatchICPPlane_LDADD		= $(LDADD)
WlzMatchICPPlane_LDFLAGS	= $(AM_LFLAGS)

WlzObjCompareSpecial_01_SOURCES	= WlzObjCompareSpecial_01.c
WlzObjCompareSpecial_01_LDADD	= $(LDADD)
WlzObjCompareSpecial_01_LDFLAGS	= $(AM_LFLAGS)

WlzOccupancy_SOURCES		= WlzOccupancy.c
WlzOccupancy_LDADD		= $(LDADD)
WlzOccupancy_LDFLAGS		= $(AM_LFLAGS)

WlzPatchObjRegister_SOURCES	= WlzPatchObjRegister.c
WlzPatchObjRegister_LDADD	= $(LDADD)
WlzPatchObjRegister_LDFLAGS	= $(AM_LFLAGS)

WlzPosterProcessForStackupKeep3StraitLines_SOURCES = \
			           WlzPosterProcessForStackupKeep3StraitLines.c
WlzPosterProcessForStackupKeep3StraitLines_LDADD = $(LDADD)
WlzPosterProcessForStackupKeep3StraitLines_LDFLAGS = $(AM_LFLAGS)

WlzProductTransform_SOURCES	= WlzProductTransform.c
WlzProductTransform_LDADD	= $(LDADD)
WlzProductTransform_LDFLAGS	= $(AM_LFLAGS)

WlzRawToWlz_SOURCES		= WlzRawToWlz.c
WlzRawToWlz_LDADD		= $(LDADD)
WlzRawToWlz_LDFLAGS		= $(AM_LFLAGS)

WlzReconstruct_SOURCES		= WlzReconstruct.c
WlzReconstruct_LDADD		= $(LDADD)
WlzReconstruct_LDFLAGS		= $(AM_LFLAGS)

WlzSizeSelect_SOURCES		= WlzSizeSelect.c
WlzSizeSelect_LDADD		= $(LDADD)
WlzSizeSelect_LDFLAGS		= $(AM_LFLAGS)

WlzSplitImage_SOURCES		= WlzSplitImage.c
WlzSplitImage_LDADD		= $(LDADD)
WlzSplitImage_LDFLAGS		= $(AM_LFLAGS)

WlzSplitMontage_SOURCES		= WlzSplitMontage.c
WlzSplitMontage_LDADD		= $(LDADD)
WlzSplitMontage_LDFLAGS		= $(AM_LFLAGS)

WlzTensor_SOURCES		= WlzTensor.c
WlzTensor_LDADD			= $(LDADD)
WlzTensor_LDFLAGS		= $(AM_LFLAGS)

WlzTensorToVtk_SOURCES		= WlzTensorToVtk.c
WlzTensorToVtk_LDADD		= $(LDADD)
WlzTensorToVtk_LDFLAGS		= $(AM_LFLAGS)

WlzTiffStackToPatches_SOURCES	= WlzTiffStackToPatches.c
WlzTiffStackToPatches_LDADD	= $(LDADD)
WlzTiffStackToPatches_LDFLAGS	= $(AM_LFLAGS)

WlzTiffStackToShade_SOURCES	= WlzTiffStackToShade.c
WlzTiffStackToShade_LDADD	= $(LDADD)
WlzTiffStackToShade_LDFLAGS	= $(AM_LFLAGS)

WlzVtxDistance_SOURCES		= WlzVtxDistance.c
WlzVtxDistance_LDADD		= $(LDADD)
WlzVtxDistance_LDFLAGS		= $(AM_LFLAGS)

WlzVtxToObj_SOURCES		= WlzVtxToObj.c
WlzVtxToObj_LDADD		= $(LDADD)
WlzVtxToObj_LDFLAGS		= $(AM_LFLAGS)
