##
# \file         Makefile.am
# \author       Bill Hill
# \date         November 2007
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
# \ingroup	LibWlzBnd
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

lib_LTLIBRARIES 	= libWlzBnd.la

libWlzBnd_la_LDFLAGS 	= -release @PACKAGE_VERSION@

libWlzBnd_la_SOURCES 	= \
			  WlzBnd3DWarp.c \
			  WlzBndBasisFn.c \
			  WlzBndFunction.c \
			  WlzBndGreyInvert.c \
			  WlzBndJavaUtils.c

include_HEADERS 	= \
			  WlzBnd.h \
			  WlzBndProto.h \
			  WlzBndType.h
