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
          $(TARGETPATH)/WSetterConstants.class \
          $(TARGETPATH)/WSetterModel.class \
	  $(TARGETPATH)/SliderRangeModel.class \
	  $(TARGETPATH)/LimitListener.class \
	  $(TARGETPATH)/LimitEvent.class

all:	$(OBJECTS)

$(TARGETPATH)/WSetterGUI.class:	$(SOURCEPATH)/WSetterGUI.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WSetterGUI.java

$(TARGETPATH)/WSetter.class:	$(SOURCEPATH)/WSetter.java \
                                $(TARGETPATH)/WSetterGUI.class
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WSetter.java

$(TARGETPATH)/WSetterConstants.class:	$(SOURCEPATH)/WSetterConstants.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WSetterConstants.java

$(TARGETPATH)/WSetterModel.class:	$(SOURCEPATH)/WSetterModel.java
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/WSetterModel.java

$(TARGETPATH)/SliderRangeModel.class:	$(SOURCEPATH)/SliderRangeModel.java 
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/SliderRangeModel.java

$(TARGETPATH)/LimitListener.class:	$(SOURCEPATH)/LimitListener.java 
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/LimitListener.java

$(TARGETPATH)/LimitEvent.class:	$(SOURCEPATH)/LimitEvent.java 
	$(JAVAC) $(FLAGS) $(SOURCEPATH)/LimitEvent.java

clean: clobber

clobber:
	-rm -rf classes/wsetter *.jar

jar: $(OBJECTS)
	cd ./classes; $(JAR) cfm ../wsetter.jar bean.tmp wsetter
