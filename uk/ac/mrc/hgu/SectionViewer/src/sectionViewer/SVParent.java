package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.io.*;
import javax.help.*;

/**
 *   <b>Must</b> be implemented by any application that uses SectionViewer components.
 *   It is implemented by the utility class SVParent2D.
 *   Alternatively it may be implemented by the main application class
 *   as in JWlzViewer.
 *   @see sectionViewer.SVParent2D
 */
public interface SVParent {

   /* methods called from the SectionViewer object */

   /**
    *   Supplies the model for the underlying 3D Woolz object.
    *   @return the WlzObjModel for this SectionViewer.
    */
   public WlzObjModel		getOBJModel();

   /**
    *   Supplies the title for a SectionViewer's Frame / InternalFrame.
    *   The title is normally an indication of Pitch, Yaw and Roll angles.
    *   @return the SectionViewer title.
    */
   public String		getSVTitleText();

   /**
    *   Supplies the object which will build the anatomy menus.
    *   @return the AnatomyBuilder object for this SectionViewer.
    */
   public AnatomyBuilder	getAnatomyBuilder();

   /**
    *   Supplies the collection of elements that make up
    *   the AnatKey.
    *   @return the Vector of AnatomyElements for this SectionViewer.
    */
   public Vector		getAnatomyElements();

   /**
    *   Supplies the collection of currently open SectionViewers.
    *   @return the Vector of currently open SectionViewers.
    */
   public Vector		getOpenViews();

   /**
    *   Removes a SectionViewer from the collection of
    *   currently open SectionViewers.
    *   @param sv the SectionViewer to be removed.
    */
   public void			removeView(SectionViewer sv);

   /**
    *   Supplies the SectionViewer Help Broker.
    *   @return the HelpBroker for SectionViewer.
    */
   public HelpBroker		getSVHelpBroker();

   /**
    *   Tests for the existence of an AnatKey.
    *   @return true if there is an AnatKey.
    */
   public AnatKey		getAnatomyKey();

   /**
    *   Returns the SVLocks object, required for synchronization.
    *   @return _Locks.
    */
   public SVLocks getSVLocks(); 
}

