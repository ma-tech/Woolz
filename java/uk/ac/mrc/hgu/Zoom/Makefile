# JAVAHEAPSIZE = 50m

SOURCEPATH = ./src/zoom/
TARGETPATH = ./classes/zoom/
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar

SOURCES = \
          $(SOURCEPATH)LimitEvent.java \
          $(SOURCEPATH)LimitListener.java \
          $(SOURCEPATH)ZoomGUI.java \
          $(SOURCEPATH)ZoomModel.java \
          $(SOURCEPATH)Zoom.java

OBJECTS = \
          $(TARGETPATH)LimitEvent.class \
          $(TARGETPATH)LimitListener.class \
          $(TARGETPATH)ZoomGUI.class \
          $(TARGETPATH)ZoomModel.class \
          $(TARGETPATH)Zoom.class

all:		$(OBJECTS) jar

classes:	$(OBJECTS)

$(OBJECTS):	$(SOURCES)
	cd ./src/zoom; $(JAVAC) *.java
	mv $(SOURCEPATH)*.class $(TARGETPATH)

jar: 		$(OBJECTS)
	cd ./classes;$(JAR) cfm ../zoom.jar bean.tmp zoom

# Use doxygen to make the documentation
doc:
	cd src/zoom; doxygen

clean: 		clobber

clobber:
	-rm -rf *.jar META-INF dox
	cd src/zoom; rm -rf *.class doxygen.tag
	cd classes; rm -rf zoom/*.class

