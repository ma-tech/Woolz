include $(top_srcdir)/config/config

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

lib_LTLIBRARIES 	= libhguDlpList.la

libhguDlpList_la_LDFLAGS = -release @PACKAGE_VERSION@

libhguDlpList_la_SOURCES = \
			  HGUDlpList.c

include_HEADERS 	= \
			  HGUDlpList.h
