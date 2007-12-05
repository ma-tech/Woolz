include $(top_srcdir)/config/config

AM_CPPFLAGS		= -I. \
			  -I$(top_srcdir)/libWlzExtFF \
			  -I$(top_srcdir)/libhguDlpList \
			  -I$(top_srcdir)/libbibfile \
			  -I$(top_srcdir)/libWlz \
			  -I$(top_srcdir)/libAlg \
			  -I$(top_srcdir)/libAlc \
             		  $(all_includes)

LDADD 			= -L$(top_srcdir)/libWlzExtFF/.libs -lWlzExtFF \
			  -L$(top_srcdir)/libhguDlpList/.libs -lhguDlpList \
			  -L$(top_srcdir)/libbibfile/.libs -lbibfile \
			  -L$(top_srcdir)/libWlz/.libs -lWlz \
			  -L$(top_srcdir)/libAlg/.libs -lAlg \
			  -L$(top_srcdir)/libAlc/.libs -lAlc \
			  -ljpeg \
			  -ltiff \
			  -lm

bin_PROGRAMS		= \
			  WlzExtFFConstruct3D \
			  WlzExtFFContourObj \
			  WlzExtFFConvert \
			  WlzExtFFWlzToXYZ \
			  WlzSAToWlz \
			  WlzTiff2Wlz

WlzExtFFConstruct3D_SOURCES	= WlzExtFFConstruct3D.c
WlzExtFFConstruct3D_LDADD	= $(LDADD)
WlzExtFFConstruct3D_LDFLAGS	= $(AM_LFLAGS)

WlzExtFFContourObj_SOURCES	= WlzExtFFContourObj.c
WlzExtFFContourObj_LDADD	= $(LDADD)
WlzExtFFContourObj_LDFLAGS	= $(AM_LFLAGS)

WlzExtFFConvert_SOURCES		= WlzExtFFConvert.c
WlzExtFFConvert_LDADD		= $(LDADD)
WlzExtFFConvert_LDFLAGS		= $(AM_LFLAGS)

WlzExtFFWlzToXYZ_SOURCES	= WlzExtFFWlzToXYZ.c
WlzExtFFWlzToXYZ_LDADD		= $(LDADD)
WlzExtFFWlzToXYZ_LDFLAGS	= $(AM_LFLAGS)

WlzSAToWlz_SOURCES		= WlzSAToWlz.c
WlzSAToWlz_LDADD		= $(LDADD)
WlzSAToWlz_LDFLAGS		= $(AM_LFLAGS)

WlzTiff2Wlz_SOURCES		= WlzTiff2Wlz.c
WlzTiff2Wlz_LDADD		= $(LDADD)
WlzTiff2Wlz_LDFLAGS		= $(AM_LFLAGS)

