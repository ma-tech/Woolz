include ../../Makefile.conf

CLASSES		= $(shell ls uk/ac/mrc/hgu/Wlz/*.class)

all:		jar

jar:
	$(JAR) cf uk.jar $(CLASSES)

clean: 		clobber

clobber:
	-rm -rf *.jar META-INF

