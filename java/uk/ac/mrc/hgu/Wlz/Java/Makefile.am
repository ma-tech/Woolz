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
WLZTOP                  = ../../../../../../..
CLASSPATH		= $(WLZTOP)/java
JAVAROOT		= $(WLZTOP)/java

JAR			= JWlz.jar

JWlzClass		= WlzObject.class

HandcodedJ              =  \
			  WlzBase.java \
			  WlzDBox2.java \
			  WlzDBox3.java \
			  WlzDVertex2.java \
			  WlzDVertex3.java \
			  WlzException.java \
			  WlzFBox2.java \
			  WlzFBox3.java \
			  WlzFileInputStream.java \
			  WlzFileOutputStream.java \
			  WlzFileStdErrStream.java \
			  WlzFileStdInStream.java \
			  WlzFileStdOutStream.java \
			  WlzFileStream.java \
			  WlzFVertex2.java \
			  WlzFVertex3.java \
			  WlzGreyV.java \
			  WlzIBox2.java \
			  WlzIBox3.java \
			  WlzIVertex2.java \
			  WlzIVertex3.java \
			  WlzNative.java \
			  WlzObject.java \
			  WlzPixelP.java \
			  WlzPixelV.java \
			  WlzPointer.java

include ../C2Java/Binding/JWlzEnumList.am 	# defines JWlzEnumList

include ../C2Java/Binding/JWlzStructList.am	# defines JWlzStructList

JWlzdir			= .
JWlz_JAVA		= \
			  $(JWlzJ)  \
			  $(HandcodedJ) \
			  $(JWlzEnumList) \
			  $(JWlzStructList)

all:			jar

jar:			$(JWlzClass)
			cd $(JAVAROOT); \
			jar cf $(JAR) uk/ac/mrc/hgu/Wlz/*.class

$(JWlzClass):		$(JWlz_JAVA)

clean:			
			$(RM) *.java *.javapp ../*.class \
			      classJWlz.stamp $(JAVAROOT)/$(JAR)

WlzBase.java:
			$(CP) ../C2Java/Binding/Enum/*.java .
			$(CP) ../C2Java/Binding/Struct/*.java .
			$(CP) ../HandCodedJava/*.java .
			$(CP) ../HandCodedJava/*.javapp .

WlzObject.java:		WlzBase.java WlzObject.javapp
			cpp -C -I. WlzObject.javapp | \
			grep -v '^#' > WlzObject.java
endif
