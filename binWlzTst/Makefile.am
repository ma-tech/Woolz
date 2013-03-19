##
# \file         binWlzTst/Makefile.am
# \author       Bill Hill
# \date         February 2012
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2012],
# The University Court of the University of Edinburgh,
# Old College, Edinburgh, UK.
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
# \brief        Process this file with automake to produce Makefile.in
# \ingroup	BinWlzTst
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
			  WlzTstCMeshCellStats \
			  WlzTstCMeshDist \
			  WlzTstCMeshGen \
			  WlzTstCMeshTransformObj \
			  WlzTstCMeshVtxInMesh \
			  WlzTstDistC \
			  WlzTstGeomArcLength2D \
			  WlzTstGeomLineTriangleIntersect \
			  WlzTstGeomLSqOPlane \
			  WlzTstGeomRectFromWideLine \
			  WlzTstGeomTetraAffineSolve \
			  WlzTstGeomTriangleAffineSolve \
			  WlzTstItrSpiral \
			  WlzTstLBTDomain \
			  WlzTstObjectCache \
			  WlzTstRegCCor \
			  WlzTstThreshold \
			  WlzTstTiledValues \
			  WlzTstVxInSimplex \
			  WlzTstGeomVtxOnLineSegment


WlzTstCMeshCellStats_SOURCES		= WlzTstCMeshCellStats.c
WlzTstCMeshCellStats_LDADD		= $(LDADD)
WlzTstCMeshCellStats_LDFLAGS		= $(AM_LFLAGS)

WlzTstCMeshDist_SOURCES			= WlzTstCMeshDist.c
WlzTstCMeshDist_LDADD			= $(LDADD)
WlzTstCMeshDist_LDFLAGS			= $(AM_LFLAGS)

WlzTstCMeshGen_SOURCES			= WlzTstCMeshGen.c
WlzTstCMeshGen_LDADD			= $(LDADD)
WlzTstCMeshGen_LDFLAGS			= $(AM_LFLAGS)

WlzTstCMeshTransformObj_SOURCES		= WlzTstCMeshTransformObj.c
WlzTstCMeshTransformObj_LDADD		= $(LDADD)
WlzTstCMeshTransformObj_LDFLAGS		= $(AM_LFLAGS)

WlzTstCMeshVtxInMesh_SOURCES		= WlzTstCMeshVtxInMesh.c
WlzTstCMeshVtxInMesh_LDADD		= $(LDADD)
WlzTstCMeshVtxInMesh_LDFLAGS		= $(AM_LFLAGS)

WlzTstDistC_SOURCES			= WlzTstDistC.c
WlzTstDistC_LDADD			= $(LDADD)
WlzTstDistC_LDFLAGS			= $(AM_LFLAGS)

WlzTstGeomArcLength2D_SOURCES		= WlzTstGeomArcLength2D.c
WlzTstGeomArcLength2D_LDADD		= $(LDADD)
WlzTstGeomArcLength2D_LDFLAGS		= $(AM_LFLAGS)

WlzTstGeomLineTriangleIntersect_SOURCES	= WlzTstGeomLineTriangleIntersect.c
WlzTstGeomLineTriangleIntersect_LDADD	= $(LDADD)
WlzTstGeomLineTriangleIntersect_LDFLAGS	= $(AM_LFLAGS)

WlzTstGeomLSqOPlane_SOURCES		= WlzTstGeomLSqOPlane.c
WlzTstGeomLSqOPlane_LDADD		= $(LDADD)
WlzTstGeomLSqOPlane_LDFLAGS		= $(AM_LFLAGS)

WlzTstGeomRectFromWideLine_SOURCES	= WlzTstGeomRectFromWideLine.c
WlzTstGeomRectFromWideLine_LDADD	= $(LDADD)
WlzTstGeomRectFromWideLine_LDFLAGS	= $(AM_LFLAGS)

WlzTstGeomTetraAffineSolve_SOURCES	= WlzTstGeomTetraAffineSolve.c
WlzTstGeomTetraAffineSolve_LDADD	= $(LDADD)
WlzTstGeomTetraAffineSolve_LDFLAGS	= $(AM_LFLAGS)

WlzTstGeomTriangleAffineSolve_SOURCES	= WlzTstGeomTriangleAffineSolve.c
WlzTstGeomTriangleAffineSolve_LDADD	= $(LDADD)
WlzTstGeomTriangleAffineSolve_LDFLAGS	= $(AM_LFLAGS)

WlzTstItrSpiral_SOURCES			= WlzTstItrSpiral.c
WlzTstItrSpiral_LDADD			= $(LDADD)
WlzTstItrSpiral_LDFLAGS			= $(AM_LFLAGS)

WlzTstLBTDomain_SOURCES			= WlzTstLBTDomain.c
WlzTstLBTDomain_LDADD			= $(LDADD)
WlzTstLBTDomain_LDFLAGS			= $(AM_LFLAGS)

WlzTstObjectCache_SOURCES		= WlzTstObjectCache.c
WlzTstObjectCache_LDADD			= $(LDADD)
WlzTstObjectCache_LDFLAGS		= $(AM_LFLAGS)

WlzTstRegCCor_SOURCES			= WlzTstRegCCor.c
WlzTstRegCCor_LDADD			= $(LDADD)
WlzTstRegCCor_LDFLAGS			= $(AM_LFLAGS)

WlzTstThreshold_SOURCES			= WlzTstThreshold.c
WlzTstThreshold_LDADD			= $(LDADD)
WlzTstThreshold_LDFLAGS			= $(AM_LFLAGS)

WlzTstTiledValues_SOURCES		= WlzTstTiledValues.c
WlzTstTiledValues_LDADD			= $(LDADD)
WlzTstTiledValues_LDFLAGS		= $(AM_LFLAGS)

WlzTstVxInSimplex_SOURCES		= WlzTstVxInSimplex.c
WlzTstVxInSimplex_LDADD			= $(LDADD)
WlzTstVxInSimplex_LDFLAGS		= $(AM_LFLAGS)

WlzTstGeomVtxOnLineSegment_SOURCES	= WlzTstGeomVtxOnLineSegment.c
WlzTstGeomVtxOnLineSegment_LDADD	= $(LDADD)
WlzTstGeomVtxOnLineSegment_LDFLAGS	= $(AM_LFLAGS)


