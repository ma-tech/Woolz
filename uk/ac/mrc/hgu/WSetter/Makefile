# JAVAHEAPSIZE = 50m

# wsetter classes

SOURCEPATH = ./src/wsetter/
TARGETPATH = ./classes/wsetter/
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar
JAVACFLAGS = -d classes

OBJECTS = $(TARGETPATH)/WSetterGUI.class \
          $(TARGETPATH)/WSetter.class \
          $(TARGETPATH)/WSetterConstants.class \
          $(TARGETPATH)/WSetterModel.class \
	  $(TARGETPATH)/SliderRangeModel.class \
	  $(TARGETPATH)/LimitListener.class \
	  $(TARGETPATH)/LimitEvent.class

all:	$(OBJECTS)

rule to make classes from .java files
$(TARGETPATH)%.class:   $(SOURCEPATH)%.java
	$(JAVAC) $(JAVACFLAGS) $<

clean: clobber

clobber:
	-rm -rf classes/wsetter/*.class *.jar META-INF

jar: $(OBJECTS)
	cd ./classes; $(JAR) cfm ../wsetter.jar bean.tmp wsetter
