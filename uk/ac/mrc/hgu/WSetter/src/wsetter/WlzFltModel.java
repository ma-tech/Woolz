//#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzFltModel.java
* \author       Nick Burton
* \date         April 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Model for a control for setting float values
* \todo         -
* \bug          None known.
*/

package wsetter;

import wsetter.*;
import javax.sound.sampled.*;
import javax.swing.event.*;

public class WlzFltModel extends FloatControl {

//----------------------------------------------------
   static class WlzType extends javax.sound.sampled.FloatControl.Type {
  	public WlzType(String str) {
	  super(str);
	}
   }

//----------------------------------------------------
   // specific types
   static WlzType DIST = new WlzType("Distance");
   static WlzType PITCH = new WlzType("Pitch");
   static WlzType YAW = new WlzType("Yaw");
   static WlzType ROLL = new WlzType("Roll");
   // generic types
   static WlzType FLOAT = new WlzType("Float");
   static WlzType INT = new WlzType("Int");

//----------------------------------------------------
   public WlzFltModel(FloatControl.Type type,
		      float min,
		      float max) {

      super (type, min, max, 0.1f, 10, min, "");
   }

//----------------------------------------------------
   // keep track of all the listeners to this model
   protected EventListenerList changeListeners =
                             new EventListenerList();

  // add a listener to the register
  public void addChangeListener(ChangeListener x) {
    changeListeners.add (ChangeListener.class, x);

    // bring it up to date with current state
    x.stateChanged(new ChangeEvent(this));
  }

  // remove a listener from the register
  public void removeChangeListener(ChangeListener x) {
    changeListeners.remove (ChangeListener.class, x);
  }


//----------------------------------------------------
   // override setValue() so that an event is fired
   public void setValue(int val) {
      if(val < super.getMinimum()) val = (int) (super.getMinimum()) + 1;
      if(val > super.getMaximum()) val = (int) (super.getMaximum()) - 1;
      super.setValue((float)val);
      fireChange();
   }

   public void setValue(float val) {
      if(val < super.getMinimum()) val = super.getMinimum();
      if(val > super.getMaximum()) val = super.getMaximum();
      super.setValue(val);
      fireChange();
   }

//----------------------------------------------------
   protected void fireChange() {
   // Create the event:
   ChangeEvent ce = new ChangeEvent(this);
   // Get the listener list
   Object[] listeners =
     changeListeners.getListenerList();
   // Process the listeners last to first
   // List is in pairs, Class and instance
   for (int i
     = listeners.length-2; i >= 0; i -= 2) {
     if (listeners[i] == ChangeListener.class) {
        ChangeListener cl = (ChangeListener)listeners[i+1];
        cl.stateChanged(ce);
     }
   }
  }

//----------------------------------------------------
   /*
   public float getMin() {
      return super.getMinimum();
   }

   public float getMax() {
      return super.getMaximum();
   }
   */
//----------------------------------------------------
   protected void debug() {
      Float min = new Float(super.getMinimum());
      Float max = new Float(super.getMaximum());

      System.out.println("debug info for WlzFltModel");
      System.out.println("min = "+min.toString());
      System.out.println("max = "+max.toString());

   }

} // class WLzFltModel
