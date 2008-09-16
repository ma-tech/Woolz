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
			  WlzTstAlgCrossCorr \
			  WlzTstCMeshDist \
			  WlzTstCMeshGen \
			  WlzTstCMeshTransformObj \
			  WlzTstGeomArcLength2D \
			  WlzTstGeomRectFromWideLine \
			  WlzTstGeomTetraAffineSolve \
			  WlzTstGeomTriangleAffineSolve \
			  WlzTstItrSpiral \
			  WlzTstLBTDomain \
			  WlzTstRegCCor \
			  WlzTstVxInSimplex \
			  WlzTstGeomVtxOnLineSegment


WlzTstAlgCrossCorr_LDADD		= $(LDADD)
WlzTstAlgCrossCorr_LDFLAGS		= $(AM_LFLAGS)

WlzTstCMeshGen_SOURCES			= WlzTstCMeshGen.c
WlzTstCMeshGen_LDADD			= $(LDADD)
WlzTstCMeshGen_LDFLAGS			= $(AM_LFLAGS)

WlzTstCMeshTransformObj_SOURCES		= WlzTstCMeshTransformObj.c
WlzTstCMeshTransformObj_LDADD		= $(LDADD)
WlzTstCMeshTransformObj_LDFLAGS		= $(AM_LFLAGS)

WlzTstGeomArcLength2D_SOURCES		= WlzTstGeomArcLength2D.c
WlzTstGeomArcLength2D_LDADD		= $(LDADD)
WlzTstGeomArcLength2D_LDFLAGS		= $(AM_LFLAGS)

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

WlzTstRegCCor_SOURCES			= WlzTstRegCCor.c
WlzTstRegCCor_LDADD			= $(LDADD)
WlzTstRegCCor_LDFLAGS			= $(AM_LFLAGS)

WlzTstVxInSimplex_SOURCES		= WlzTstVxInSimplex.c
WlzTstVxInSimplex_LDADD			= $(LDADD)
WlzTstVxInSimplex_LDFLAGS		= $(AM_LFLAGS)

WlzTstGeomVtxOnLineSegment_SOURCES	= WlzTstGeomVtxOnLineSegment.c
WlzTstGeomVtxOnLineSegment_LDADD	= $(LDADD)
WlzTstGeomVtxOnLineSegment_LDFLAGS	= $(AM_LFLAGS)


