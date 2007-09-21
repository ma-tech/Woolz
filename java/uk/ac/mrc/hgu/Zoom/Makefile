# JAVAHEAPSIZE = 50m

SOURCEPATH = ./src/zoom/
TARGETPATH1 = ./classes/
TARGETPATH2 = ./classes/zoom/
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar
#JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1) -verbose
#JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1) -g
JAVACFLAGS = -classpath $(CLASSPATH) -sourcepath $(SOURCEPATH) -d $(TARGETPATH1)
JAVADOC = /net/mahost/export/system0/java/j2sdk1.4.1_01/bin/javadoc
JAVADOCFLAGS = -sourcepath ./src -d docs

all:		zoom jar

zoom:
	$(JAVAC) $(JAVACFLAGS) $(SOURCEPATH)*.java

jar:
	-cp -rf images $(TARGETPATH1)
	cd $(TARGETPATH1); $(JAR) cf ../zoom.jar zoom images

# Use doxygen to make the documentation
doc:
	cd $(SOURCEPATH); doxygen

# Use javadoc to make the documentation
javadocs:
	$(JAVADOC) $(JAVADOCFLAGS) zoom

clean: 		clobber

clobber:
	-rm -rf *.jar META-INF dox docs
	cd $(SOURCEPATH); rm -rf *.class doxygen.tag
	cd $(TARGETPATH2); rm -rf *.class

