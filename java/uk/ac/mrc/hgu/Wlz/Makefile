#########################################################################
# Project:      Woolz
# Title:        Makefile
# Date:         August 1999
# Author:       Bill Hill
# Copyright:	1999 Medical Research Council, UK.
#		All rights reserved.
# Address:	MRC Human Genetics Unit,
#		Western General Hospital,
#		Edinburgh, EH4 2XU, UK.
# Purpose:      Makefile for the Java/Woolz binding.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../../../../../../Makefile.conf

HANDCODEDJAVA		= $(shell cd HandCodedJava; ls *.java)
CLASSPATH		= ../../../../..
JAVACFLAGS		= -classpath $(CLASSPATH)
JAVA_PACKAGE		= uk.ac.mrc.hgu.Wlz
JARFILE			= Wlz.jar

all:			binding archive classes jar

clean:
			$(RM) *.java *.class
			cd C; $(MAKE) clean
			cd C2Java; $(MAKE) clean

install:		install_archive install_jar

binding:
			cd C2Java; $(MAKE) JAVA_PACKAGE=$(JAVA_PACKAGE) binding

ifeq ($(UNIXTYPE), CYGWIN) 
archive:
			$(CP) C2Java/Binding/WlzObjectC2J.c C
			$(CP) C2Java/Binding/WlzObjectC2J.h C
			$(CP) HandCodedC/*.c HandCodedC/*.h C
			cd C; $(MAKE) dll
else
archive:
			$(CP) C2Java/Binding/WlzObjectC2J.c C
			$(CP) C2Java/Binding/WlzObjectC2J.h C
			$(CP) HandCodedC/*.c HandCodedC/*.h C
			cd C; $(MAKE) soarchive
endif							

ifeq ($(UNIXTYPE), CYGWIN) 
install_archive:
			cd C; $(MAKE) install_dll
else
install_archive:
			cd C; $(MAKE) install_soarchive
endif

install_jar:		$(JARFILE)
			cd ../../../../..; \
			$(INSTALL) -o $(OWNER) -m 644 $(JARFILE) $(HGU_LIB_DIR)
			

$(JARFILE):		WlzObject.class
			cd ../../../../..; \
			jar -cf $(JARFILE) uk/ac/mrc/hgu/Wlz/*.class

classes:		WlzObject.class

WlzObject.class:	WlzObject.java
			$(JAVAC)  $(JAVACFLAGS) *.java

WlzObject.java:
			$(CP) C2Java/Binding/Enum/*.java \
			      C2Java/Binding/Struct/*.java . ; \
			for J in $(HANDCODEDJAVA) ; \
			do \
			  $(CPP) -C -I. HandCodedJava/$$J | \
			  grep -v '^#' >$$J ; \
			done
jar:
			cd ../../../../..;gmake clean jar
