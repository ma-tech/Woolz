# JAVAHEAPSIZE = 50m

SOURCEPATH = ./src/zoom
TARGETPATH = ./classes/zoom
#JAVAC = /opt/java/bin/javac
JAVAC = /net/mahost/export/system0/java/j2sdk1.4.1_01/bin/javac
JAR = /opt/java/bin/jar
JAVACFLAGS = -d classes

OBJECTS = $(TARGETPATH)/ZoomGUI.class \
          $(TARGETPATH)/Zoom.class \
          $(TARGETPATH)/ZoomModel.class \
	  $(TARGETPATH)/LimitListener.class \
	  $(TARGETPATH)/LimitEvent.class

all:	$(OBJECTS)

# rule to make classes from .java files
$(TARGETPATH)%.class:   $(SOURCEPATH)%.java
	$(JAVAC) $(JAVACFLAGS) $<

clean: clobber

clobber:
	-rm -rf classes/zoom/*.class classes/*.class *.jar META-INF

jar: $(OBJECTS)
	cd ./classes;cp -r images zoom;$(JAR) cfm ../zoom.jar bean.tmp zoom
