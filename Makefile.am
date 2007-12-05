## Makefile.am -- Process this file with automake to produce Makefile.in

include $(top_srcdir)/config/config

SUBDIRS =	\
		libAlc \
		libAlg \
		libWlz \
		binWlz \
		binWlzTst
if BUILD_EXTFF
  SUBDIRS +=	\
		libbibfile \
		libhguDlpList \
		libWlzExtFF \
		libReconstruct \
		binWlzApp \
		binWlzExtFF
endif

doc:
		doxygen Doxyfile_Core

install-doc:	doc
		cp -r html_Core ${docdir}/WlzCore

