package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

import uk.ac.mrc.hgu.Wlz.*;

/**
 *   Corresponds to a single row of the Anatomy Key.
 *   (ie a KeyElement object).
 *   @see sectionViewer.AnatKey
 */
public class AnatomyElement {

   
   /**   The 3D anatomy component.  */
   private WlzObject _obj;
   /**   The full name of the anatomy component.  */
   private String _descr;

   /**   A unique index number for the anatomy component.  */
   private int _indx;

   /**   True if the anatomy component will be displayed.  */
   private boolean _visible;
   /**   True if the anatomy component has been removed from the AnatKey.  */
   //private boolean _removed;

   /**
    *   Creates an AnatomyElement with the given anatomy component (obj), 
    *   and name (descr).
    *   @param obj A 3D Woolz object representing an anatomy component.
    *   @param descr A String giving the full name of the anatomy component.
    */
   public AnatomyElement(WlzObject obj,
                         String descr,
			 int indx) {
      _obj = obj;
      _descr = descr;
      _indx = indx;
      _visible = true;
      //_removed = false;
   }

   /**
    *   Returns the anatomy component.
    *   @return the 3D Woolz object representing the anatomy component.
    */
   public WlzObject getObj() {
      return _obj;
   }

   /**
    *   Returns the unique index number.
    *   @return _indx
    */
   public int getIndx() {
      return _indx;
   }

   /**
    *   Returns true if the anatomy component is visible.
    *   @return true if anatomy component is visible.
    *   Visible here means that the anatomy component will be displayed
    *   if it is intersected by the 2D section.
    */
   public boolean isVisible() {
      return _visible;
   }

   /**
    *   Toggles the visibility of the anatomy component.
    *   @param state true if anatomy component is visible.
    */
   public void setElementVisible(boolean state) {
      _visible = state;
   }

   /**
    *   Returns the removed or replaced status of the anatomy component.
    *   @return true if anatomy component has been removed.
    */
/*
   public boolean isRemoved() {
      return _removed;
   }
*/

   /**
    *   Toggles the removed or replaced status of the anatomy component.
    *   @param state true if anatomy component is to be removed.
    */
/*
   public void setRemoved(boolean state) {
      _removed = state;
   }
*/

   /**
    *   Returns the full name of the anatomy component.
    *   @return the full name of the anatomy component.
    */
   protected String getDescription() {
      return _descr;
   }

   /**
    *   Sets the full name of the anatomy component.
    *   @param descr the full name of the anatomy component.
    */
   protected void setDescription(String descr) {
      _descr = descr;
   }


}

 
