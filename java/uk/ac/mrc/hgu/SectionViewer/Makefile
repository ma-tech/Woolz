# Makefile for SectionViewer

include ../../../../../../../Makefile.conf

CLASSPATH = .:..:../../../../..:../WSetter/classes:../Zoom/classes

ifeq ($(UNIXTYPE), CYGWIN) 
CLASSPATH = '.;..;../../../../..;../WSetter/classes;../Zoom/classes'
endif
SOURCEPATH = ./src/sectionViewer/
TARGETPATH1 = ./classes/
TARGETPATH2 = ./classes/sectionViewer/
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar
#JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1) -verbose
#JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1) -g
JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1)
JAVADOC = /net/mahost/export/system0/java/j2sdk1.4.1_01/bin/javadoc
JAVADOCFLAGS = -sourcepath ./src -d docs

ifeq ($(UNIXTYPE), DARWIN)
	JAVAC = /usr/bin/javac
	JAR = /usr/bin/jar
endif

all:	sectionViewer jar

sectionViewer:
	$(JAVAC) $(JAVACFLAGS) $(SOURCEPATH)*.java

jar:	
	-cp -rf images $(TARGETPATH1)
	cd $(TARGETPATH1); $(JAR) cf ../sectionViewer.jar sectionViewer images

# Use doxygen to make the documentation
dox:
	cd $(SOURCEPATH); doxygen

# Use javadoc to make the documentation
javadocs:
	$(JAVADOC) $(JAVADOCFLAGS) sectionViewer

clean: clobber

clobber:
	-rm -rf *.jar META-INF dox docs
	cd $(SOURCEPATH); rm -rf *.class doxygen.tag
	cd $(TARGETPATH2); rm -rf *.class

