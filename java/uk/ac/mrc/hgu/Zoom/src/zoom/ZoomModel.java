package zoom;

import zoom.*;
import java.util.*;
import javax.swing.event.*;

/**
 * Model class for Zoom bean.
 * <br>Uses the <b>Model View Controller</b> paradigm.
 * @author Nick Burton
 * @see Zoom
 */
public class ZoomModel {

private int _modelMin = 20;
private int _modelMax = 1000;
private int _modelVal = 100;
private int _modelInc = 50;

//----------------------------------------------------
   public ZoomModel() {
   }

//----------------------------------------------------
   public int getValue() {
      return _modelVal;
   }
//----------------------------------------------------
   public void setValue(int val) {
      _modelVal = val;
      fireChange();
   }
//----------------------------------------------------
   public void incValue(boolean inc) {
      if(inc) {
	 _modelVal += _modelInc;
	 if(_modelVal > _modelMax) _modelVal = _modelMax;
      } else {
	 _modelVal -= _modelInc;
	 if(_modelVal < _modelMin) _modelVal = _modelMin;
      }
      fireChange();
   }
//----------------------------------------------------
   public int getMin() {
      return _modelMin;
   }
//----------------------------------------------------
   public void setMin(int val) {
      _modelMin = val;
   }
//----------------------------------------------------
   public int getMax() {
      return _modelMax;
   }
//----------------------------------------------------
   public void setMax(int val) {
      _modelMax = val;
   }
//----------------------------------------------------
   public int getInc() {
      return _modelInc;
   }
//----------------------------------------------------
   public void setInc(int val) {
      _modelInc = val;
   }
//====================================================
  /**
   * Fires an event if the model value changes.
   * @param     void
   * @return    void
   */
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
   // keep track of all the listeners to this model
   protected EventListenerList changeListeners =
                             new EventListenerList();

  // add a listener to the register
  /**
   * Adds a ChangeListener for events fired from this bean.
   * @param     ChangeListener x
   * @return    void
   */
  public void addChangeListener(ChangeListener x) {
    changeListeners.add (ChangeListener.class, x);

    // bring it up to date with current state
    x.stateChanged(new ChangeEvent(this));
  }

  // remove a listener from the register
  /**
   * Removes a ChangeListener for events fired from this bean.
   * @param     ChangeListener x
   * @return    void
   */
  public void removeChangeListener(ChangeListener x) {
    changeListeners.remove (ChangeListener.class, x);
  }

//----------------------------------------------------
  protected synchronized void fireLimitChange() {
     // Create the event:
     LimitEvent le = new LimitEvent(this);
     // Iterate through the listener list
     Iterator listeners = _limitListeners.iterator();
     while(listeners.hasNext()) {
	((LimitListener)listeners.next()).limitChanged(le);
     }
  } // fireLimitChange

//----------------------------------------------------
  // keep track of all the LimitListeners to this model
  protected List _limitListeners = new ArrayList();

  public synchronized void addLimitListener(LimitListener x) {
     _limitListeners.add (x);
  }

  public synchronized void removeLimitListener(LimitListener x) {
     _limitListeners.remove (x);
  }

} // class ZoomModel
