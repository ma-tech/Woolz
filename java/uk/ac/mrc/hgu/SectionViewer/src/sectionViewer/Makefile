# JAVAHEAPSIZE = 150m

# viewer classes

include ../../../../Makefile.conf

TARGETPATH = ../../classes/sectionViewer/
JAVAC = /opt/java/bin/javac
JAVACFLAGS = -g


all:
	$(JAVAC) $(JAVACFLAGS) *.java
	mv *.class $(TARGETPATH)

clean: clobber

clobber:
	rm -rf *.class

