# JAVAHEAPSIZE = 50m

# wsetter classes

CLASSPATH = .:./src
SOURCEPATH = ./src/wsetter
TARGETPATH = ./classes/wsetter
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar
FLAGS = -classpath $(CLASSPATH) -d classes

OBJECTS = $(TARGETPATH)/WSetter.class \
          $(TARGETPATH)/WlzFltModel.class \
	  $(TARGETPATH)/SliderRangeModel.class \
	  $(TARGETPATH)/test.class

all:	$(OBJECTS)

$(TARGETPATH)/WSetter.class:	$(SOURCEPATH)/WSetter.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WSetter.java

$(TARGETPATH)/WlzFltModel.class:	$(SOURCEPATH)/WlzFltModel.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WlzFltModel.java

$(TARGETPATH)/SliderRangeModel.class:	$(SOURCEPATH)/SliderRangeModel.java 
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/SliderRangeModel.java

$(TARGETPATH)/test.class:	$(SOURCEPATH)/test.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/test.java


clean: clobber

clobber:
	-rm -rf classes/wsetter *.jar

jarfile: $(OBJECTS)
	cd ./classes; $(JAR) cfm ../wsetter.jar bean.tmp wsetter
