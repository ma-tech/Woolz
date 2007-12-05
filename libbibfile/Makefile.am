include $(top_srcdir)/config/config

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

lib_LTLIBRARIES 	= libbibfile.la

libbibfile_la_LDFLAGS 	= -release @PACKAGE_VERSION@

libbibfile_la_SOURCES 	= \
			  bibFileAlloc.c \
			  bibFileIO.c \
			  bibFileParse.c

include_HEADERS 	= \
			  bibFile.h
