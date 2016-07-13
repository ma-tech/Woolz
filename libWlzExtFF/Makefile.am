##
# \file         libWlzExtFF/Makefile.am
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
# \bug          None known.
##

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libhguDlpList \
			  -I$(top_srcdir)/libbibfile \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
			  -I$(includedir) \
			  $(INCL_EXTFF) \
             		  $(all_includes)

lib_LTLIBRARIES 	= libWlzExtFF.la

libWlzExtFF_la_LDFLAGS 	= -release @PACKAGE_VERSION@

libWlzExtFF_la_SOURCES 	= \
			  WlzExtFFAm.c \
			  WlzExtFFAnl.c \
			  WlzExtFFBibUtils.c \
			  WlzExtFFBmp.c \
			  WlzExtFF.c \
			  WlzExtFFDen.c \
			  WlzExtFFEMT.c \
			  WlzExtFFGif.c \
			  WlzExtFFIcs.c \
			  WlzExtFFIPL.c \
			  WlzExtFFJpeg.c \
			  WlzExtFFMesh.c \
			  WlzExtFFNifti.c \
			  WlzExtFFNodeEle.c \
			  WlzExtFFNrrd.c \
			  WlzExtFFObj.c \
			  WlzExtFFPic.c \
			  WlzExtFFPly2.c \
			  WlzExtFFPnm.c \
			  WlzExtFFPvl.c \
			  WlzExtFFSlc.c \
			  WlzExtFFSMesh.c \
			  WlzExtFFStack.c \
			  WlzExtFFStl.c \
			  WlzExtFFTiff.c \
			  WlzExtFFTxt.c \
			  WlzExtFFVff.c \
			  WlzExtFFVMesh.c \
			  WlzExtFFVtk.c

#			  WlzExtFFRaw.c 

include_HEADERS 	= \
			  WlzExtFF.h \
			  WlzExtFFType.h \
			  WlzExtFFProto.h
