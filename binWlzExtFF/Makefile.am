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
			  -I$(top_srcdir)/libWlzExtFF \
			  -I$(top_srcdir)/libhguDlpList \
			  -I$(top_srcdir)/libbibfile \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

LDADD 			= -L$(top_srcdir)/libWlzExtFF/.libs -lWlzExtFF \
			  -L$(top_srcdir)/libhguDlpList/.libs -lhguDlpList \
			  -L$(top_srcdir)/libbibfile/.libs -lbibfile \
			  -L$(top_srcdir)/libWlz/.libs -lWlz \
			  -L$(top_srcdir)/libAlg/.libs -lAlg \
			  -L$(top_srcdir)/libAlc/.libs -lAlc \
			  ${LIBS_EXTFF} ${LIBS}

bin_PROGRAMS		= \
			  WlzExtFFConstruct3D \
			  WlzExtFFContourObj \
			  WlzExtFFConvert \
			  WlzExtFFWlzToXYZ \
			  WlzSAToWlz \
			  WlzTiff2Wlz

WlzExtFFConstruct3D_SOURCES	= WlzExtFFConstruct3D.c
WlzExtFFConstruct3D_LDADD	= $(LDADD)
WlzExtFFConstruct3D_LDFLAGS	= $(AM_LFLAGS)

WlzExtFFContourObj_SOURCES	= WlzExtFFContourObj.c
WlzExtFFContourObj_LDADD	= $(LDADD)
WlzExtFFContourObj_LDFLAGS	= $(AM_LFLAGS)

WlzExtFFConvert_SOURCES		= WlzExtFFConvert.c
WlzExtFFConvert_LDADD		= $(LDADD)
WlzExtFFConvert_LDFLAGS		= $(AM_LFLAGS)

WlzExtFFWlzToXYZ_SOURCES	= WlzExtFFWlzToXYZ.c
WlzExtFFWlzToXYZ_LDADD		= $(LDADD)
WlzExtFFWlzToXYZ_LDFLAGS	= $(AM_LFLAGS)

WlzSAToWlz_SOURCES		= WlzSAToWlz.c
WlzSAToWlz_LDADD		= $(LDADD)
WlzSAToWlz_LDFLAGS		= $(AM_LFLAGS)

WlzTiff2Wlz_SOURCES		= WlzTiff2Wlz.c
WlzTiff2Wlz_LDADD		= $(LDADD)
WlzTiff2Wlz_LDFLAGS		= $(AM_LFLAGS)

