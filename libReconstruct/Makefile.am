include $(top_srcdir)/config/config

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libWlzExtFF \
			  -I$(top_srcdir)/libhguDlpList \
			  -I$(top_srcdir)/libbibfile \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

lib_LTLIBRARIES 	= libReconstruct.la

libReconstruct_la_LDFLAGS  = -release @PACKAGE_VERSION@

libReconstruct_la_SOURCES  = \
			  ReconstructAuto.c \
			  ReconstructConstruct3D.c \
			  ReconstructCrossCor.c \
			  ReconstructDebug.c \
			  ReconstructExplode3D.c \
			  ReconstructFileIO.c \
			  ReconstructMisc.c \
			  ReconstructPreProc.c \
			  ReconstructRegister.c \
			  ReconstructRotMatch.c \
			  ReconstructSection.c \
			  ReconstructTranMatch.c

include_HEADERS 	= \
			  Reconstruct.h \
			  ReconstructProto.h \
			  ReconstructType.h
