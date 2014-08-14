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
JAVASRCS		= \
			  WlzAffineTransformLSq.java \
			  WlzArrayStats.java \
			  WlzArrayViewer.java \
			  WlzBndGreyInvert.java \
			  WlzBoundingBoxIntersection.java \
			  WlzBoundingBox.java \
			  WlzConstruct3D.java \
			  WlzDitherViewer.java \
			  WlzExtFFFacts.java \
			  WlzFacts.java \
			  WlzGetBackground.java \
			  WlzGreyStats.java \
			  WlzLabelTest.java \
			  WlzTest3DVtx.java \
			  WlzUnion.java \
			  WlzViewer.java

CLASSPATH		= .:../../../../../..

JAVAFLAGS		= -classpath $(CLASSPATH)

JAVACLS			= $(JAVASRCS:.java=.class)

.java.class:
			$(JAVAC) $(JAVAFLAGS) $<

all:			$(JAVACLS)

clean:
			$(RM) *.class

demo:			$(JAVACLS)
			@echo "Simple demos of using Woolz from Java."
			@echo "Look at the corresponding Java source file to"
			@echo "see the code!"
			@echo
			@echo "All tests will require the JWlz dynamic libary"
			@echo "to be in the library search path."
			@echo
			@echo "WlzBoundingBox simply prints the bounding box"
			@echo "of the given woolz object:"
			$(JAVA) $(JAVAFLAGS) WlzBoundingBox t0.wlz
			@echo 
			@echo "WlzBoundingBoxIntersection prints the"
			@echo "intersection of the bounding box of the"
			@echo "given woolz object with a section plane:"
			$(JAVA) $(JAVAFLAGS) WlzBoundingBoxIntersection \
					     -d 20.0 -p 90.0 t0.wlz
			@echo 
			@echo "WlzDitherViewer throws up a window with the"
			@echo "second object dithered onto the first (the"
			@echo "window needs to be resized):"
			$(JAVA) $(JAVAFLAGS) WlzDitherViewer t0.wlz t1.wlz
			@echo 
			@echo "WlzFacts prints basic facts about the woolz"
			@echo "object:"
			$(JAVA) $(JAVAFLAGS) WlzFacts t0.wlz
			@echo 
			@echo "WlzGetBackground prints the objects background"
			@echo "value"
			$(JAVA) $(JAVAFLAGS) WlzGetBackground t0.wlz
			@echo 
			@echo "WlzGreyStats prints simple grey values stats"
			@echo "computed from the given object:"
			$(JAVA) $(JAVAFLAGS) WlzGreyStats t0.wlz
			@echo 
			@echo "WlzArrayStats calculates simple statistics on"
			@echo "the values in a 2D or 3D array created from"
			@echo "a Woolz object."
			$(JAVA) $(JAVAFLAGS) WlzArrayStats t0.wlz
			@echo 
			@echo "WlzArrayViewer throws up a window with the"
			@echo "the given object. (the window needs to be"
			@echo "resized):"
			$(JAVA) $(JAVAFLAGS) WlzArrayViewer t0.wlz
			@echo 
			@echo "WlzUnion computes the union of objects given"
			@echo "on the command line:"
			$(JAVA) $(JAVAFLAGS) WlzUnion -o out.wlz t0.wlz t1.wlz
			@echo 
			@echo "WlzViewer throws up a window with the"
			@echo "the given object, single left mouse click"
			@echo "shows grey value (the window needs to be"
			@echo "resized):"
			$(JAVA) $(JAVAFLAGS) WlzViewer t0.wlz
endif
