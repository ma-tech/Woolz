package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

import uk.ac.mrc.hgu.Wlz.*;

/**
 *   Corresponds to a single row of the Anatomy Key
 *   @see sectionViewer.AnatKey
 */
public class AnatomyElement {

   
   /**   The 3D anatomy component.  */
   private WlzObject _obj;
   /**   The full name of the anatomy component.  */
   private String _descr;

   /**   The position (row) of the anatomy component in the AnatKey.  */
   private int _indx;
   /**   True if the anatomy component will be displayed.  */
   private boolean _visible;
   /**   True if the anatomy component has been removed from the AnatKey.  */
   private boolean _removed;

   /**
    *   Creates an AnatomyElement with the given anatomy component (obj), 
    *   name (descr) and position in the AnatKey (indx).
    *   @param obj A 3D Woolz object representing an anatomy component.
    *   @param descr A String giving the full name of the anatomy component.
    *   @param indx An int giving the position (row) of the anatomy component
    *   in the AnatKey. Possible values are 0 to 5.
    */
   public AnatomyElement(WlzObject obj, String descr, int indx) {
      _obj = obj;
      _descr = descr;
      _indx = indx;
      _visible = true;
      _removed = false;
   }

   /**
    *   Returns the anatomy component.
    *   @param
    *   @return the 3D Woolz object representing the anatomy component.
    */
   protected WlzObject getObj() {
      return _obj;
   }

   /**
    *   Returns true if the anatomy component is visible.
    *   @param
    *   @return true if anatomy component is visible.
    *   Visible here means that the anatomy component will be displayed
    *   if it is intersected by the 2D section.
    */
   protected boolean isVisible() {
      return _visible;
   }

   /**
    *   Toggles the visibility of the anatomy component.
    *   @param state true if anatomy component is visible.
    */
   public void setVisible(boolean state) {
      _visible = state;
   }

   /**
    *   Returns the removed or replaced status of the anatomy component.
    *   @param
    *   @return true if anatomy component has been removed.
    */
   public boolean isRemoved() {
      return _removed;
   }

   /**
    *   Toggles the removed or replaced status of the anatomy component.
    *   @param state true if anatomy component is to be removed.
    */
   public void setRemoved(boolean state) {
      _removed = state;
   }

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

   /**
    *   Returns the position (row) of the anatomy component
    *   in the AnatKey.
    *   @return the position (row) of the anatomy component.
    */
   protected int getIndex() {
      return _indx;
   }

   /**
    *   Sets the position (row) of the anatomy component
    *   in the AnatKey.
    *   @param indx the position of the anatomy component.
    */
   protected void setIndex(int indx) {
      _indx = indx;
   }

   /**
    *   Returns the position (row) of the next available
    *   row in the AnatKey.
    *   If an anatomy component has been removed from the AnatKey
    *   this row is immediately available.
    *   If no anatomy components have been deleted the next available 
    *   row is the next empty one (or row 0 if all 6 rows are occupied).
    *   @param arr an array of AnatomyElements from which the next available 
    *   space is deduced.
    *   @return the next available row in the AnatKey..
    */
   public static int getNextIndex(AnatomyElement[] arr) {
      // return the lowest available empty space in the array

      AnatomyElement el = null;
      int ret = 0;

      for(int i=0; i<AnatKey._nrows; i++) {
         el = arr[i];
         if((el == null) || (el.isRemoved())) {
	    ret = i;
	    break;
	 }
      }
      return ret;
   }
}

 
