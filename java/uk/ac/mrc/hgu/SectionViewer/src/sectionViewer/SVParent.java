package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.io.*;
import javax.help.*;

public interface SVParent {

   /*
      methods that need to be implemented by any class
      that uses SectionViewer objects
    */

   /* methods called from the SectionViewer object */
   public WlzObjModel		getOBJModel();

   public String		getSVTitleText();

   public AnatomyBuilder	getAnatomyBuilder();

   public AnatomyElement[]	getAnatomyArr();

   public Vector		getOpenViews();

   public void			removeView(SectionViewer sv);

   public HelpBroker		getSVHelpBroker();

}

