# JAVAHEAPSIZE = 50m

CLASSPATH = .:./src
SOURCEPATH = ./src/zoom
TARGETPATH = ./classes/zoom
#JAVAC = /opt/java/bin/javac
JAVAC = /net/mahost/export/system0/java/j2sdk1.4.1_01/bin/javac
JAR = /opt/java/bin/jar
FLAGS = -classpath $(CLASSPATH) -d classes

OBJECTS = $(TARGETPATH)/ZoomGUI.class \
          $(TARGETPATH)/Zoom.class \
          $(TARGETPATH)/ZoomModel.class \
	  $(TARGETPATH)/LimitListener.class \
	  $(TARGETPATH)/LimitEvent.class

all:	$(OBJECTS)

$(TARGETPATH)/ZoomGUI.class:	$(SOURCEPATH)/ZoomGUI.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/ZoomGUI.java

$(TARGETPATH)/Zoom.class:	$(SOURCEPATH)/Zoom.java \
                                $(TARGETPATH)/ZoomGUI.class
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/Zoom.java

$(TARGETPATH)/ZoomModel.class:	$(SOURCEPATH)/ZoomModel.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/ZoomModel.java

$(TARGETPATH)/LimitListener.class:	$(SOURCEPATH)/LimitListener.java 
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/LimitListener.java

$(TARGETPATH)/LimitEvent.class:	$(SOURCEPATH)/LimitEvent.java 
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/LimitEvent.java

clean: clobber

clobber:
	-rm -rf classes/zoom classes/*.class *.jar

jar: $(OBJECTS)
	cd ./classes;cp -r images zoom;$(JAR) cfm ../zoom.jar bean.tmp zoom
