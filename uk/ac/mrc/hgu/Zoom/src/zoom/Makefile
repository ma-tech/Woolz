# JAVAHEAPSIZE = 50m

include ../../../../../../../../../Makefile.conf	# up 9 levels

#JAVAC = /opt/java/bin/javac
JAVAC = /net/mahost/export/system0/java/j2sdk1.4.1_01/bin/javac
JAR = /opt/java/bin/jar
JAVACFLAGS = 

OBJECTS = \
	  LimitEvent.class \
	  LimitListener.class \
          ZoomGUI.class \
          ZoomModel.class \
          Zoom.class

all:	$(OBJECTS)

clean: clobber

clobber:
	-rm -rf *.class
