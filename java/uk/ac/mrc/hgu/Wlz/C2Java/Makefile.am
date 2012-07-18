##
# \file         Makefile.am
# \author       Bill Hill
# \date         July 2012
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2012],
# The University Court of the University of Edinburgh,
# Old College, Edinburgh, UK.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief        Process this file with automake to produce Makefile.in
##

if BUILD_JAVA
WLZTOP			= ../../../../../../../

WLZINCDIRS		= \
			-I$(WLZTOP)/libAlc \
			-I$(WLZTOP)/libAlg \
			-I$(WLZTOP)/libWlz \
			-I$(WLZTOP)/libWlzExtFF \
			-I$(WLZTOP)/libbibfile \
			-I$(WLZTOP)/libWlzBnd

WLZINCLIDES		= \
			$(WLZTOP)libAlc/Alc.h \
			$(WLZTOP)libAlc/AlcType.h \
			$(WLZTOP)libAlc/AlcProto.h \
			$(WLZTOP)libAlg/Alg.h \
			$(WLZTOP)libAlg/AlgType.h \
			$(WLZTOP)libAlg/AlgProto.h \
			$(WLZTOP)libWlz/Wlz.h \
			$(WLZTOP)libWlz/WlzError.h \
			$(WLZTOP)libWlz/WlzType.h \
			$(WLZTOP)libWlz/WlzDebug.h \
			$(WLZTOP)libWlz/WlzProto.h \
			$(WLZTOP)libWlz/WlzMacro.h \
			$(WLZTOP)libbibfile/bibFile.h \
			$(WLZTOP)libWlzExtFF/WlzExtFF.h \
			$(WLZTOP)libWlzExtFF/WlzExtFFType.h \
			$(WLZTOP)libWlzExtFF/WlzExtFFProto.h \
			$(WLZTOP)libWlzBnd/WlzBnd.h \
			$(WLZTOP)libWlzBnd/WlzBndType.h \
			$(WLZTOP)libWlzBnd/WlzBndProto.h

all:		 	binding

clean: 
			$(RM) *.c \
			      *.class \
			      timestamp \
			      Binding/Enum/* \
			      Binding/Struct/* \
			      Binding/WlzObjectC2J.javapp \
			      Binding/WlzObjectC2J.c  \
			      Binding/WlzObjectC2J.h

clobber:		clean
			$(RM)  *.java

binding:		parser timestamp

parser:			WlzC2Java.class


timestamp:		WlzC2Java.c
			$(MKDIR_P) Binding/Enum
			$(MKDIR_P) Binding/Struct
			@JAVA@ WlzC2Java -package uk.ac.mrc.hgu.Wlz \
				WlzC2Java.c
			touch timestamp


WlzC2Java.java:		WlzC2Java.jj
			@JAVACC@ WlzC2Java.jj

WlzC2Java.class:	WlzC2Java.java
			@JAVAC@ WlzC2Java.java

WlzC2Java.c:		WlzC2Java.cpp $(WLZINCLIDES)
			$(CC) -E $(WLZINCDIRS) -DWLZ_EXT_BIND \
			      WlzC2Java.cpp > WlzC2Java.c
endif
