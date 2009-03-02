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
             		  $(all_includes)

lib_LTLIBRARIES 	= libAlc.la

libAlc_la_LDFLAGS 	= -release @PACKAGE_VERSION@

libAlc_la_SOURCES 	= \
			  AlcAlloc.c \
			  AlcArray.c \
			  AlcBlockStack.c \
			  AlcCPQueue.c \
			  AlcDLPList.c \
			  AlcFreeStack.c \
			  AlcHashTable.c \
			  AlcHeap.c \
			  AlcKDTree.c \
			  AlcString.c \
			  AlcVector.c

include_HEADERS 	= \
			  Alc.h \
			  AlcProto.h \
			  AlcTemplates.h \
			  AlcType.h
