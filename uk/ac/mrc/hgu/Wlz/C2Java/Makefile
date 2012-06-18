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
# Purpose:      Sub-makefile for the Java/Woolz binding.
# $Revision$
# $Id$
# Maintenance:	Log changes below, with most recent at top of list.
#########################################################################

# Include the configuration makefile (modify as required).
include			../../../../../../../../Makefile.conf

# Name of archive (library) to be built (modify as required).
ARCHIVE			=
SOARCHIVE		=

# Names of executables to be built (modify as required).
EXECUTABLES		=

# List of all 'C' source files (modify as required).
CSOURCES		=

# List of all 'C pre-processor' source files (modify as required).
CPPSOURCES		= WlzC2Java.cpp

# List of all 'Java' source files (modify as required).
JAVASOURCES		=

# List of all 'JavaCC' source files (modify as required).
JAVACCSOURCES		= WlzC2Java.jj

# Java package (modify as required).
JAVA_PACKAGE		= uk.ac.mrc.hgu.Wlz

# List of all header files that are available outside of either this archive
# or these executables (modify as required).
INCLUDES_PUB		=

# List of all header files that are required by this archive or these
# executables (modify as required).
INCLUDES_PRV		=

# List of all troff sources for manual pages, these should have the appropriate
# extension (eg cat.1) (modify as required).
MANPAGES		=

# List of all directories to be searched before the standard include
# and lib directories (modify as required).
DEPENDDIRS		= \
			  ../../../../../../../libAlc \
			  ../../../../../../../libAlg \
			  ../../../../../../../libWlz \
			  ../../../../../../../libWlzExtFF \
			  ../../../../../../../libbibfile \
			  ../../../../../../../libWlzBnd

# List of all header file directories to be searched (modify as required).
INCDIRS			= \
			  $(DEPENDDIRS) \
			  $(HGU_INC_DIR)

# List of library search paths (modify as required).
LIBDIRS			= \
			  $(DEPENDDIRS) \
			  $(SYSLIB) \
			  $(HGU_LIB_DIR)

# List of libraries to link (modify as required).

# List of all sub directories that contain makefiles. The subsystems will be
# built in the order in which they are given here (modify as required).
SUBSYSTEMS              =

# List of all header files, ie the union of the private and public header
# files (should not need modifying).
INCLUDES_ALL		= $(sort $(INCLUDES_PUB) $(INCLUDES_PRV))

# List of all object files required for archive (should not need modifying).
OBJECTS			= $(CSOURCES:%.c=%.o)

# Basic flags for controlling compilation (modify as required).
DEFINES         	= -DWLZ_EXT_BIND $(THREADS) \
			  -D$(UNIXTYPE) $(UNIXFLAGS) -DHGUHOME=\"$(HGUHOME)\"
CPPFLAGS		= $(INCDIRS:%=-I%) $(DEFINES)

# Target which builds everything
all:			binding

# Target which builds the Java binding
binding:		parser parse-timestamp

# Target for the Java parser
parser:			WlzC2Java.class

parse-timestamp:	WlzC2Java.c
			$(JAVA) WlzC2Java -package $(JAVA_PACKAGE) WlzC2Java.c
			touch parse-timestamp

WlzC2Java.c:
			$(CPP) $(CPPFLAGS) \
			$(CPPSOURCES) > WlzC2Java.c

WlzC2Java.class:	WlzC2Java.java
			$(JDK_HOME)/bin/javac WlzC2Java.java

WlzC2Java.java:		WlzC2Java.jj


# Target to clean up everything
clobber:		clean

# Target to clean up after making
clean:			
			-$(RM) *.c *.class *.java \
			parse-timestamp Binding/Enum/* \
			Binding/Struct/* \
			Binding/WlzObjectC2J.*
