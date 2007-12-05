include $(top_srcdir)/config/config

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
			  WlzTstCMeshGen \
			  WlzTstCMeshTransformObj \
			  WlzTstGeomTetraAffineSolve \
			  WlzTstGeomTriangleAffineSolve \
			  WlzTstItrSpiral \
			  WlzTstLBTDomain \
			  WlzTstRegCCor \
			  WlzTstVxInSimplex

Wlz3DWarpMQ_SOURCES			= WlzTstAlgCrossCorr.c
WlzTstAlgCrossCorr_LDADD		= $(LDADD)
WlzTstAlgCrossCorr_LDFLAGS		= $(AM_LFLAGS)

WlzTstCMeshGen_SOURCES			= WlzTstCMeshGen.c
WlzTstCMeshGen_LDADD			= $(LDADD)
WlzTstCMeshGen_LDFLAGS			= $(AM_LFLAGS)

WlzTstCMeshTransformObj_SOURCES		= WlzTstCMeshTransformObj.c
WlzTstCMeshTransformObj_LDADD		= $(LDADD)
WlzTstCMeshTransformObj_LDFLAGS		= $(AM_LFLAGS)

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

