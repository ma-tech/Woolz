# JAVAHEAPSIZE = 150m

include ../../../../../../../Makefile.conf

CLASSPATH = .:..:../../../../../../..:../../../WSetter/classes:../../../Zoom/classes
SOURCEPATH = ./src/sectionViewer/
TARGETPATH = ./classes/sectionViewer/
JAVAC = /opt/java/bin/javac
JAR = /opt/java/bin/jar
JAVACFLAGS = -classpath $(CLASSPATH)

SOURCES = \
	  $(SOURCEPATH)SelectableMenu.java \
	  $(SOURCEPATH)ScrollableTextField.java \
	  $(SOURCEPATH)PointEntryGUI.java \
	  $(SOURCEPATH)PointEntry.java \
	  $(SOURCEPATH)AnatKey.java \
	  $(SOURCEPATH)AnatKeyGUI.java \
	  $(SOURCEPATH)AnatomyElement.java \
	  $(SOURCEPATH)AnatomyBuilder.java \
	  $(SOURCEPATH)ParseXML.java \
	  $(SOURCEPATH)WlzImgView.java \
	  $(SOURCEPATH)WlzObjModel.java \
	  $(SOURCEPATH)ViewStructModel.java \
	  $(SOURCEPATH)SVParent.java \
	  $(SOURCEPATH)SectionViewerGUI.java \
	  $(SOURCEPATH)SectionViewer.java \
	  $(SOURCEPATH)SVFileFilter.java \
	  $(SOURCEPATH)SVUtils.java \
	  $(SOURCEPATH)SVStatus.java

OBJECTS = \
	  $(TARGETPATH)SelectableMenu.class \
	  $(TARGETPATH)ScrollableTextField.class \
	  $(TARGETPATH)PointEntryGUI.class \
	  $(TARGETPATH)PointEntry.class \
	  $(TARGETPATH)AnatKey.class \
	  $(TARGETPATH)AnatKeyGUI.class \
	  $(TARGETPATH)AnatomyElement.class \
	  $(TARGETPATH)AnatomyBuilder.class \
	  $(TARGETPATH)ParseXML.class \
	  $(TARGETPATH)WlzImgView.class \
	  $(TARGETPATH)WlzObjModel.class \
	  $(TARGETPATH)ViewStructModel.class \
	  $(TARGETPATH)SVParent.class \
	  $(TARGETPATH)SectionViewerGUI.class \
	  $(TARGETPATH)SectionViewer.class \
	  $(TARGETPATH)SVFileFilter.class \
	  $(TARGETPATH)SVUtils.class \
	  $(TARGETPATH)SVStatus.class



all:	classes jar

classes:	$(OBJECTS)

$(OBJECTS):	$(SOURCES)
	cd $(SOURCEPATH); $(JAVAC) $(JAVACFLAGS) *.java
	mv $(SOURCEPATH)*.class $(TARGETPATH)

jar:	classes
	cd classes; $(JAR) cfm ../sectionViewer.jar bean.tmp sectionViewer

# Use doxygen to make the documentation
doc:
	cd $(SOURCEPATH); doxygen

clean: clobber

clobber:
	-rm -rf *.jar META-INF dox
	cd $(SOURCEPATH); rm -rf *.class doxygen.tag
	cd $(TARGETPATH); rm -rf *.class

