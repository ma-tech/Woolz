# Makefile for WSetter

SOURCEPATH = ./src/wsetter/
TARGETPATH1 = ./classes/
TARGETPATH2 = ./classes/wsetter/
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar
#JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1) -verbose
#JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1) -g
JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1)
JAVADOC = /net/mahost/export/system0/java/j2sdk1.4.1_01/bin/javadoc
JAVADOCFLAGS = -sourcepath ./src -d docs

all:		wsetter jar

wsetter:
	$(JAVAC) $(JAVACFLAGS) $(SOURCEPATH)*.java

jar:
	cd $(TARGETPATH1); $(JAR) cf ../wsetter.jar wsetter

# Use doxygen to make the documentation
doc:
	cd $(SOURCEPATH); doxygen

# Use javadoc to make the documentation
javadocs:
	$(JAVADOC) $(JAVADOCFLAGS) wsetter

clean: 		clobber

clobber:
	-rm -rf *.jar META-INF dox docs
	cd $(SOURCEPATH); rm -rf *.class doxygen.tag
	cd $(TARGETPATH2); rm -rf *.class

