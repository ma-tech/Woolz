package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

import uk.ac.mrc.hgu.Wlz.*;

public class AnatomyElement {

   private WlzObject _obj;
   private String _descr;

   private int _indx;
   private boolean _visible;
   private boolean _removed;

   public AnatomyElement(WlzObject obj, String descr, int indx) {
      _obj = obj;
      _descr = descr;
      _indx = indx;
      _visible = true;
      _removed = false;
   }

   protected WlzObject getObj() {
      return _obj;
   }

   protected boolean isVisible() {
      return _visible;
   }

   public void setVisible(boolean state) {
      _visible = state;
   }

   public boolean isRemoved() {
      return _removed;
   }

   public void setRemoved(boolean state) {
      _removed = state;
   }

   protected String getDescription() {
      return _descr;
   }

   protected void setDescription(String descr) {
      _descr = descr;
   }

   protected int getIndex() {
      return _indx;
   }

   protected void setIndex(int indx) {
      _indx = indx;
   }

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

 
