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
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

lib_LTLIBRARIES 	= libAlg.la

libAlg_la_LDFLAGS 	= -release @PACKAGE_VERSION@

libAlg_la_SOURCES 	= \
			  AlgAutoCorr.c \
			  AlgBits.c \
			  AlgComplexUtils.c \
			  AlgConvolve.c \
			  AlgCrossCorr.c \
			  AlgDebug.c \
			  AlgDPSearch.c \
			  AlgFourier.c \
			  AlgGamma.c \
			  AlgHeapSort.c \
			  AlgLinearFit.c \
			  AlgMatrixCG.c \
			  AlgMatrixGauss.c \
			  AlgMatrixLU.c \
			  AlgMatrixMath.c \
			  AlgMatrixRSEigen.c \
			  AlgMatrixRSTDiag.c \
			  AlgMatrixSV.c \
			  AlgMatrixTDiagQLI.c \
			  AlgMixture.c \
			  AlgPolyLSQ.c \
			  AlgQSort.c \
			  AlgRand.c \
			  AlgRandZig.c \
			  AlgRange.c \
			  AlgRank.c \
			  AlgShuffle.c \
			  AlgVectorMath.c

include_HEADERS 	= \
			  Alg.h \
			  AlgProto.h \
			  AlgType.h
