include $(top_srcdir)/config/config

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

lib_LTLIBRARIES 	= libWlz.la

libWlz_la_LDFLAGS 	= -release @PACKAGE_VERSION@

libWlz_la_SOURCES 	= \
			  WlzBnd3DWarp.c \
			  WlzBndBasisFn.c \
			  WlzBndFunction.c \
			  WlzBndGreyInvert.c \
			  WlzBndJavaUtils.c

include_HEADERS 	= \
			  WlzBnd.h \
			  WlzBndProto.h \
			  WlzBndType.h
