# JAVAHEAPSIZE = 50m

SOURCEPATH = ./src/wsetter/
TARGETPATH = ./classes/wsetter/
#JAVAC = /opt/java/bin/javac
JAVAC = /net/mahost/export/system0/java/j2sdk1.4.1_01/bin/javac
JAR = /opt/java/bin/jar

SOURCES = \
          $(SOURCEPATH)LimitEvent.java \
          $(SOURCEPATH)LimitListener.java \
          $(SOURCEPATH)WSetterConstants.java \
          $(SOURCEPATH)WSetterGUI.java \
          $(SOURCEPATH)WSetterModel.java \
          $(SOURCEPATH)SliderRangeModel.java \
          $(SOURCEPATH)WSetter.java

OBJECTS = \
          $(TARGETPATH)LimitEvent.class \
          $(TARGETPATH)LimitListener.class \
          $(TARGETPATH)WSetterConstants.class \
          $(TARGETPATH)WSetterGUI.class \
          $(TARGETPATH)WSetterModel.class \
          $(TARGETPATH)SliderRangeModel.class \
          $(TARGETPATH)WSetter.class

all:		$(OBJECTS) jar

classes:	$(OBJECTS)

$(OBJECTS):	$(SOURCES)
	cd ./src/wsetter; gmake all
	mv $(SOURCEPATH)*.class $(TARGETPATH)

jar: 		$(OBJECTS)
	cd ./classes;$(JAR) cfm ../wsetter.jar bean.tmp wsetter

# Use doxygen to make the documentation
doc:
	cd src/wsetter; doxygen

clean: 		clobber

clobber:
	-rm -rf *.jar META-INF dox
	cd src/wsetter; rm -rf *.class doxygen.tag
	cd classes; rm -rf wsetter/*.class

