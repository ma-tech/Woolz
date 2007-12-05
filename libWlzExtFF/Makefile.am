include $(top_srcdir)/config/config

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libhguDlpList \
			  -I$(top_srcdir)/libbibfile \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
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
			  WlzExtFFGif.c \
			  WlzExtFFIcs.c \
			  WlzExtFFIPL.c \
			  WlzExtFFJpeg.c \
			  WlzExtFFPic.c \
			  WlzExtFFPnm.c \
			  WlzExtFFSlc.c \
			  WlzExtFFStack.c \
			  WlzExtFFTiff.c \
			  WlzExtFFVff.c \
			  WlzExtFFVtk.c

#			  WlzExtFFRaw.c 

include_HEADERS 	= \
			  WlzExtFF.h
