##
# \file         Makefile.am
# \author       Bill Hill
# \date         July 2012
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
##

if BUILD_JAVA
AM_CPPFLAGS		= \
			  -I. \
			  -I$(top_srcdir)/libReconstruct \
			  -I$(top_srcdir)/libWlzExtFF \
			  -I$(top_srcdir)/libhguDlpList \
			  -I$(top_srcdir)/libbibfile \
			  -I$(top_srcdir)/libWlzBnd \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
			  $(all_includes)

libJWlz_la_LIBADD	= -L$(top_srcdir)/libWlzExtFF/.libs -lWlzExtFF \
			  -L$(top_srcdir)/libhguDlpList/.libs -lhguDlpList \
			  -L$(top_srcdir)/libbibfile/.libs -lbibfile \
			  -L$(top_srcdir)/libWlzBnd/.libs -lWlzBnd \
			  -L$(top_srcdir)/libWlz/.libs -lWlz \
			  -L$(top_srcdir)/libAlg/.libs -lAlg \
			  -L$(top_srcdir)/libAlc/.libs -lAlc \
			  -L$(libdir) \
			  ${LIBS_EXTFF} \
			  ${LIBS}

JWlzdir			= .

HandcodedC		= \
			  WlzFileStream.c \
			  WlzJavaArray1D.c \
			  WlzJavaArray2D.c \
			  WlzJavaArray3D.c \
			  WlzJavaUtils.c \
			  WlzJavaValue.c \
			  WlzNative.c

HandcodedH		= \
			  WlzFileStream.h \
			  WlzJavaArray1D.h \
			  WlzJavaArray2D.h \
			  WlzJavaArray3D.h \
			  WlzJava.h \
			  WlzJavaUtils.h \
			  WlzJavaValue.h

JWlzGenC		= WlzObjectC2J.c

JWlzGenH		= WlzObjectC2J.h

lib_LTLIBRARIES		= libJWlz.la

libJWlz_la_LDFLAGS	= -release @PACKAGE_VERSION@

libJWlz_la_SOURCES	= $(HandcodedC) $(JWlzGenC)

include_HEADERS		= $(HandcodedH) $(JWlzGenH)

clean:			
			$(RM) *.lo *.o *.c *.h $(lib_LTLIBRARIES)
			$(RMRF) .libs

WlzFileStream.c:
			$(CP) ../HandCodedC/*.c .
			$(CP) ../HandCodedC/*.h .
			$(CP) ../C2Java/Binding/WlzObjectC2J.c .
			$(CP) ../C2Java/Binding/WlzObjectC2J.h .
endif
