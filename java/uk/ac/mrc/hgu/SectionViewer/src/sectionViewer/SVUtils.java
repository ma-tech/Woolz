package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.io.*;
import javax.help.*;

interface SVUtils {

   /*
   methods that need to be implemented by any class
   that uses SectionViewer objects
   */

   /* these methods are called from the SectionViewer object */
  WlzObjModel		getOBJModel();

  String		getTitleText();

  AnatomyBuilder	getAnatomyBuilder();

  AnatomyElement[]	getAnatomyArr();

  Vector		getOpenViews();

  HelpBroker		getHelpBroker();

   /* these methods are called from the controlling object */
  void			openGreyLevel(File imgFile);

  void			closeGreyLevel();
}

