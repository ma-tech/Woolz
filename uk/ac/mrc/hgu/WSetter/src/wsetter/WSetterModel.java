package wsetter;

import wsetter.*;
import java.util.*;
import javax.swing.event.*;

/**
 * Model class for WSetter bean.
 * <br>Uses the <b>Model View Controller</b> paradigm.
 * @author Nick Burton
 * @see WSetter
 * @see SliderRangeModel
 */
public class WSetterModel implements WSetterConstants {

private int _type;

/*
private int minInt;
private int maxInt;
private int intVal;
private float minFlt;
private float maxFlt;
private float fltVal;
private double minDbl;
private double maxDbl;
private double dblVal;
*/

private double modelMin = 0.0;
private double modelMax = 1000.0;
private double modelInitVal = 0.0;

private Vector minVec = null;
private Vector maxVec = null;
private Vector valVec = null;
//----------------------------------------------------
  /**
   * Constructors
   */
   public WSetterModel(int type) {
      _type = type;
      minVec = new Vector(1,0);
      maxVec = new Vector(1,0);
      valVec = new Vector(1,0);
      switch(_type) {
         case INTEGER:
	    minVec.insertElementAt(new Integer((int)modelMin), 0);
	    maxVec.insertElementAt(new Integer((int)modelMax), 0);
	    valVec.insertElementAt(new Integer((int)modelInitVal), 0);
	    break;
         case FLOAT:
	    minVec.insertElementAt(new Float((float)modelMin), 0);
	    maxVec.insertElementAt(new Float((float)modelMax), 0);
	    valVec.insertElementAt(new Float((float)modelInitVal), 0);
	    break;
         case DOUBLE:
	    minVec.insertElementAt(new Double(modelMin), 0);
	    maxVec.insertElementAt(new Double(modelMax), 0);
	    valVec.insertElementAt(new Double(modelInitVal), 0);
	    break;
	 default:
	    break;
      }
   }

//----------------------------------------------------
  /**
   * @return
   */
   public Vector getValue() {
      return valVec;
   }
//----------------------------------------------------
  /**
   * @return    void
   */
   public void setValue(double val) {
      switch(_type) {
         case INTEGER:
	    valVec.setElementAt(new Integer((int)val), 0);
	    fireChange();
	    break;
         case FLOAT:
	    valVec.setElementAt(new Float((float)val), 0);
	    fireChange();
	    break;
         case DOUBLE:
	    valVec.setElementAt(new Double(val), 0);
	    fireChange();
	    break;
	 default:
	    break;
      }
   }
//----------------------------------------------------
   public Vector getMin() {
      return minVec;
   }
//----------------------------------------------------
   public void setMin(double val) {
      switch(_type) {
         case INTEGER:
	    minVec.setElementAt(new Integer((int)val), 0);
	    fireLimitChange();
	    break;
         case FLOAT:
	    minVec.setElementAt(new Float((float)val), 0);
	    fireLimitChange();
	    break;
         case DOUBLE:
	    minVec.setElementAt(new Double(val), 0);
	    fireLimitChange();
	    break;
	 default:
	    break;
      }
   }
//----------------------------------------------------
   public Vector getMax() {
      return maxVec;
   }
//----------------------------------------------------
   public void setMax(double val) {
      switch(_type) {
         case INTEGER:
	    maxVec.setElementAt(new Integer((int)val), 0);
	    fireLimitChange();
	    break;
         case FLOAT:
	    maxVec.setElementAt(new Float((float)val), 0);
	    fireLimitChange();
	    break;
         case DOUBLE:
	    maxVec.setElementAt(new Double(val), 0);
	    fireLimitChange();
	    break;
	 default:
	    break;
      }
   }
//----------------------------------------------------
   public int getType() {
      return _type;
   }
//====================================================
  /**
   * Fires an event if the slider value changes.
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

} // class WSetterModel
