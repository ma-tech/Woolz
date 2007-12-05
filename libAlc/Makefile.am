include $(top_srcdir)/config/config

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
			  AlcKDTree.c \
			  AlcString.c \
			  AlcVector.c

include_HEADERS 	= \
			  Alc.h \
			  AlcProto.h \
			  AlcTemplates.h \
			  AlcType.h
