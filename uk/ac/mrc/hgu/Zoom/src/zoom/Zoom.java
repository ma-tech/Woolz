package zoom;

import zoom.*;

import javax.swing.*;
import java.awt.*;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

/**
 * View / Controller class for Zoom.
 * @author Nick Burton
 * @see ZoomModel
 */
public class Zoom extends ZoomGUI {
//                     implements Serializable {

   private boolean _enabled;
//-------------------------------------------------------------
   public ZoomModel _mod1 = null;

   public buttonToModelAdaptor_A BToM_A = null;
   public buttonToModelAdaptor_B BToM_B = null;
   public textToModelAdaptor TToM = null;
   public modelToTextAdaptor MToT = null;
   public modelToThisAdaptor MToThis = null;

   // model properties
   private static int maxTextlen = 6;
//-------------------------------------------------------------
  /**
   * Creates a Zoom with event firing initially disabled.
   */
  public Zoom() {
    _enabled = false;
    try {
      initZoom();
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }
//-------------------------------------------------------------
  /**
   * Sets the Zoom's ability to fire events.
   * @param	boolean state; if true the Zoom can fire events
   * @return	void
   */
  public void setEnabled(boolean state) {
     _enabled = state;
  }
//-------------------------------------------------------------
  /**
   * Queries the Zoom's ability to fire events.
   * @param	void
   * @return	boolean	true if the Zoom can fire events
   */
  public boolean isEnabled() {
     return _enabled;
  }
//-------------------------------------------------------------
  /**
   * Initialises the Zoom's model and GUI.
   * @param	void
   * @return	void
   */
  private void initZoom() throws Exception {

    // the main Zoom model
    _mod1 = new ZoomModel();
    setModel();

  } // initZoom()
//-------------------------------------------------------------
  /**
   * Hooks up the event adaptors after removing any existing ones.
   * @param	void
   * @return	void
   */
  protected void setModel() {

     //-----------------------------------------
     // hook the adapters in place
     // (after removing any old ones)
     //-----------------------------------------

     // CONTROL ADAPTORS
     if(BToM_A != null) {
	_inButton.removeMouseListener(BToM_A);
     }
     BToM_A = new buttonToModelAdaptor_A();
     _inButton.addMouseListener(BToM_A);

     if(BToM_B != null) {
	_outButton.removeMouseListener(BToM_B);
     }
     BToM_B = new buttonToModelAdaptor_B();
     _outButton.addMouseListener(BToM_B);

     if(TToM != null) {
	textf.removeActionListener(TToM);
     }
     TToM = new textToModelAdaptor(textf, _mod1);
     //TToM = new textToModelAdaptor();
     textf.addActionListener(TToM);

//...............................
     // VIEW ADAPTORS
     if(MToT != null) {
	_mod1.removeChangeListener(MToT);
     }
     MToT = new modelToTextAdaptor(_mod1, textf);
    // MToT = new modelToTextAdaptor();
     _mod1.addChangeListener(MToT);

//...............................
     // MODEL ADAPTORS
     if(MToThis != null) {
	_mod1.removeChangeListener(MToThis);
     }
     MToThis = new modelToThisAdaptor();
     _mod1.addChangeListener(MToThis);

  } // setModel

//===========================================
// accessor methods for the Zoom properties etc
//===========================================
   public int getValue() {
      return _mod1.getValue();
   }
   public int getMin() {
      return _mod1.getMin();
   }
   public int getMax() {
      return _mod1.getMax();
   }
//.....................................
   public int getLabelWidth() {
      return _labelWidth;
   }
   public int getTextWidth() {
      return _textWidth;
   }
   public Color getBgc() {
      return _bgc;
   }
//.....................................
   public void setValue(int val) {
      _mod1.setValue(val);
   }
   public void setMin(int val) {
      _mod1.setMin(val);
   }
   public void setMax(int val) {
      _mod1.setMax(val);
   }
   //----------------------------------------------------
   public void setInc(int val) {
     _mod1.setInc(val);
   }
//.....................................
   public void setLabelWidth(int w) {
      super.setLabelWidth(w);
   }
   public void setTextWidth(int w) {
      super.setTextWidth(w);
   }
   public void setBgc(Color col) {
      super.setBgc(col);
   }
   public void setZoomLabel(String labl) {
      super.setZoomLabel(labl);
   }

//===========================================
// inner classes
//===========================================
// adapters for MVC model

// CONTROL ADAPTORS
//---------------------------------------
/**
 * Increments the Zoom's model when the + button is pressed
 */
  public class buttonToModelAdaptor_A
    implements MouseListener {

    public buttonToModelAdaptor_A() {
    }

    public void mouseClicked(MouseEvent e) {
    }
    public void mouseEntered(MouseEvent e) {
    }
    public void mouseExited(MouseEvent e) {
    }
    public void mousePressed(MouseEvent e) {
    }
    public void mouseReleased(MouseEvent e) {
       _mod1.incValue(true);
    }
  } // class buttonToModelAdaptor_A

//---------------------------------------
/**
 * Decrements the Zoom's model when the - button is pressed
 */
  public class buttonToModelAdaptor_B
    implements MouseListener {

    public buttonToModelAdaptor_B() {
    }

    public void mouseClicked(MouseEvent e) {
    }
    public void mouseEntered(MouseEvent e) {
    }
    public void mouseExited(MouseEvent e) {
    }
    public void mousePressed(MouseEvent e) {
    }
    public void mouseReleased(MouseEvent e) {
       _mod1.incValue(false);
    }

  } // class buttonToModelAdaptor_B

//---------------------------------------
  /**
   * Modifies the Zoom's model when a numeric value has been entered in the Zoom's text field.
   */
  public class textToModelAdaptor implements ActionListener {
     ZoomModel model;
     JTextField control;
     double VAL;
     String textstr = "";
     String msg = "the zoom control's text field expects a number\n such as 100 or 100.5";

     /**
      * Creates a textToModelAdaptor with the specified model and text field.
      * @param JTextField cntrl; the Zoom's text field
      * @param ZoomModel mdl; the Zoom's model
      */
     public textToModelAdaptor(JTextField cntrl, ZoomModel mdl) {
	model = mdl;
	control = cntrl;
     }

     /**
      * 
      */
     public textToModelAdaptor() {
     }

     /**
      * @param	ActionEvent e
      * @return	void
      */
     public void actionPerformed(ActionEvent e) {
	try {
	   textstr = control.getText();
	   VAL = Double.parseDouble(textstr);
	   model.setValue((int)VAL);
	}
	catch(NullPointerException npe) {
	   System.out.println(
		 "null pointer exception in Zoom textToModelAdaptor");
	}
	catch(NumberFormatException npe) {
	   control.setText("");
	   JOptionPane.showMessageDialog(null,
		 msg,
		 "alert",
		 JOptionPane.ERROR_MESSAGE);
	}

     }
  } // class textToModelAdaptor

//---------------------------------------
// VIEW ADAPTORS
//---------------------------------------
  /**
   * Updates the text field view of the Zoom when the Zoom's model has changed.
   */
  public class modelToTextAdaptor implements ChangeListener {
     ZoomModel model;
     JTextField view;
     String valstr;

     /**
      * Creates a modelToTextAdaptor with the specified model and text field.
      * @param ZoomModel mdl; the Zoom's model
      * @param JTextField vw; the Zoom's text field
      */
     public modelToTextAdaptor(ZoomModel mdl, JTextField vw) {
	view = vw;
	model = mdl;
     }

     /**
      * 
      */
     public modelToTextAdaptor() {
     }

     /**
      * Changes the text field when the model changes.
      * @param	ChangeEvent e
      * @return	void
      */
     public void stateChanged(ChangeEvent e) {
	valstr = Integer.toString(_mod1.getValue());

	if(valstr.length() > maxTextlen) {
	   textf.setText(valstr.substring(0,maxTextlen));
	} else {
	   textf.setText(valstr);
	}
     }
  } // class modelToTextAdaptor


//---------------------------------------
  // MODEL ADAPTORS
//---------------------------------------
  // allows this Zoom to fire events
  /**
   * Fires an event when the Zoom's model changes.
   */
  public class modelToThisAdaptor implements ChangeListener {

     /**
      * 
      */
     public modelToThisAdaptor() {
     }

     /**
      * Fires an event from the Zoom when the model changes.
      * @param	ChangeEvent e
      * @return	void
      */
     public void stateChanged(ChangeEvent e) {
	fireChange();
     }
  } // class modelToThisAdaptor

//===========================================
// fire events and manage listeners
//===========================================
  /**
   * List of all changeListeners registered with this Zoom.
   */
  protected EventListenerList changeListeners = new EventListenerList();

  /**
   * Adds a ChangeListener for events fired from this Zoom.
   * @param	ChangeListener x
   * @return	void
   */
  public void addChangeListener(ChangeListener x) {
     changeListeners.add (ChangeListener.class, x);

     // bring it up to date with current state
     /**
      * 
      */
     x.stateChanged(new ChangeEvent(this));
  }

  /**
   * Removes a ChangeListener for events fired from this Zoom.
   * @param	ChangeListener x
   * @return	void
   */
  public void removeChangeListener(ChangeListener x) {
     changeListeners.remove (ChangeListener.class, x);
  }

  /**
   * Fires an event from the Zoom.
   * @param	void
   * @return	void
   */
  protected void fireChange() {
     if(_enabled == true) {
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
  } // fireChange


} // class Zoom
