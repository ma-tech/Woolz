# JAVAHEAPSIZE = 50m

# wsetter classes

CLASSPATH = .:./src
SOURCEPATH = ./src/wsetter
TARGETPATH = ./classes/wsetter
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar
FLAGS = -classpath $(CLASSPATH) -d classes

OBJECTS = $(TARGETPATH)/WSetterGUI.class \
          $(TARGETPATH)/WSetter.class \
          $(TARGETPATH)/WlzFltModel.class \
	  $(TARGETPATH)/SliderRangeModel.class

all:	$(OBJECTS)

$(TARGETPATH)/WSetterGUI.class:	$(SOURCEPATH)/WSetterGUI.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WSetterGUI.java

$(TARGETPATH)/WSetter.class:	$(SOURCEPATH)/WSetter.java \
                                $(TARGETPATH)/WSetterGUI.class
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WSetter.java

$(TARGETPATH)/WlzFltModel.class:	$(SOURCEPATH)/WlzFltModel.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WlzFltModel.java

$(TARGETPATH)/SliderRangeModel.class:	$(SOURCEPATH)/SliderRangeModel.java 
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/SliderRangeModel.java

clean: clobber

clobber:
	-rm -rf classes/wsetter *.jar

jar: $(OBJECTS)
	cd ./classes; $(JAR) cfm ../wsetter.jar bean.tmp wsetter
