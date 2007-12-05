include $(top_srcdir)/config/config

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
			  AlgRange.c \
			  AlgRank.c \
			  AlgShuffle.c \
			  AlgVectorMath.c

include_HEADERS 	= \
			  Alg.h \
			  AlgProto.h \
			  AlgType.h
