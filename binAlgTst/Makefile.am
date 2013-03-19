##
# \file         Makefile.am
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
##

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

LDADD 			= -L$(top_srcdir)/libAlg/.libs -lAlg \
			  -L$(top_srcdir)/libAlc/.libs -lAlc \
			  -lm

bin_PROGRAMS		= \
			  AlgTstConvolve1 \
			  AlgTstCrossCorr1 \
			  AlgTstFourier \
			  AlgTstGamma1 \
			  AlgTstGrayCode \
			  AlgTstHeapSort \
			  AlgTstHilbertIndex \
			  AlgTstLinearFit1 \
			  AlgTstMatrixArithmetic1 \
			  AlgTstMatrixArithmetic2 \
			  AlgTstMatrixArithmetic3 \
			  AlgTstMatrixCGSolve1 \
			  AlgTstMatrixCGSolve2 \
			  AlgTstMatrixRSEigen1 \
			  AlgTstMatrixSolve1 \
			  AlgTstMixtureMLG1 \
			  AlgTstQSort1 \
			  AlgTstRange1 \
			  AlgTstRank1 \
			  AlgTstShuffle1 \
			  AlgTstShuffle2

AlgTstConvolve1_SOURCES			= AlgTstConvolve1.c
AlgTstConvolve1_LDADD			= $(LDADD)
AlgTstConvolve1_LDFLAGS			= $(AM_LFLAGS)

AlgTstCrossCorr1_SOURCES		= AlgTstCrossCorr1.c
AlgTstCrossCorr1_LDADD			= $(LDADD)
AlgTstCrossCorr1_LDFLAGS		= $(AM_LFLAGS)

AlgTstFourier_SOURCES			= AlgTstFourier.c
AlgTstFourier_LDADD			= $(LDADD)
AlgTstFourier_LDFLAGS			= $(AM_LFLAGS)

AlgTstGamma1_SOURCES			= AlgTstGamma1.c
AlgTstGamma1_LDADD			= $(LDADD)
AlgTstGamma1_LDFLAGS			= $(AM_LFLAGS)

AlgTstGrayCode_SOURCES			= AlgTstGrayCode.c
AlgTstGrayCode_LDADD			= $(LDADD)
AlgTstGrayCode_LDFLAGS			= $(AM_LFLAGS)

AlgTstLinearFit1_SOURCES		= AlgTstLinearFit1.c
AlgTstLinearFit1_LDADD			= $(LDADD)
AlgTstLinearFit1_LDFLAGS		= $(AM_LFLAGS)

AlgTstHeapSort_SOURCES			= AlgTstHeapSort.c
AlgTstHeapSort_LDADD			= $(LDADD)
AlgTstHeapSort_LDFLAGS			= $(AM_LFLAGS)

AlgTstHilbertIndex_SOURCES		= AlgTstHilbertIndex.c
AlgTstHilbertIndex_LDADD		= $(LDADD)
AlgTstHilbertIndex_LDFLAGS		= $(AM_LFLAGS)

AlgTstMatrixArithmetic1_SOURCES		= AlgTstMatrixArithmetic1.c
AlgTstMatrixArithmetic1_LDADD		= $(LDADD)
AlgTstMatrixArithmetic1_LDFLAGS		= $(AM_LFLAGS)

AlgTstMatrixArithmetic2_SOURCES		= AlgTstMatrixArithmetic2.c
AlgTstMatrixArithmetic2_LDADD		= $(LDADD)
AlgTstMatrixArithmetic2_LDFLAGS		= $(AM_LFLAGS)

AlgTstMatrixArithmetic3_SOURCES		= AlgTstMatrixArithmetic3.c
AlgTstMatrixArithmetic3_LDADD		= $(LDADD)
AlgTstMatrixArithmetic3_LDFLAGS		= $(AM_LFLAGS)

AlgTstMatrixCGSolve1_SOURCES		= AlgTstMatrixCGSolve1.c
AlgTstMatrixCGSolve1_LDADD		= $(LDADD)
AlgTstMatrixCGSolve1_LDFLAGS		= $(AM_LFLAGS)

AlgTstMatrixCGSolve2_SOURCES		= AlgTstMatrixCGSolve2.c
AlgTstMatrixCGSolve2_LDADD		= $(LDADD)
AlgTstMatrixCGSolve2_LDFLAGS		= $(AM_LFLAGS)

AlgTstMatrixRSEigen1_SOURCES		= AlgTstMatrixRSEigen1.c
AlgTstMatrixRSEigen1_LDADD		= $(LDADD)
AlgTstMatrixRSEigen1_LDFLAGS		= $(AM_LFLAGS)

AlgTstMatrixSolve1_SOURCES		= AlgTstMatrixSolve1.c
AlgTstMatrixSolve1_LDADD		= $(LDADD)
AlgTstMatrixSolve1_LDFLAGS		= $(AM_LFLAGS)

AlgTstMixtureMLG1_SOURCES		= AlgTstMixtureMLG1.c
AlgTstMixtureMLG1_LDADD			= $(LDADD)
AlgTstMixtureMLG1_LDFLAGS		= $(AM_LFLAGS)

AlgTstQSort1_SOURCES			= AlgTstQSort1.c
AlgTstQSort1_LDADD			= $(LDADD)
AlgTstQSort1_LDFLAGS			= $(AM_LFLAGS)

AlgTstRange1_SOURCES			= AlgTstRange1.c
AlgTstRange1_LDADD			= $(LDADD)
AlgTstRange1_LDFLAGS			= $(AM_LFLAGS)

AlgTstRank1_SOURCES			= AlgTstRank1.c
AlgTstRank1_LDADD			= $(LDADD)
AlgTstRank1_LDFLAGS			= $(AM_LFLAGS)

AlgTstShuffle1_SOURCES			= AlgTstShuffle1.c
AlgTstShuffle1_LDADD			= $(LDADD)
AlgTstShuffle1_LDFLAGS			= $(AM_LFLAGS)

AlgTstShuffle2_SOURCES			= AlgTstShuffle2.c
AlgTstShuffle2_LDADD			= $(LDADD)
AlgTstShuffle2_LDFLAGS			= $(AM_LFLAGS)
