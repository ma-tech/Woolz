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
 * View / Controller class for zoom control bean.
 * <br>Uses the <b>Model View Controller</b> paradigm.
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
   * Constructor
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
   * Allow or inhibit the component's ability to fire events
   * @param	state	true if the component can fire events
   * @return	void
   */
  public void setEnabled(boolean state) {
     _enabled = state;
  }
//-------------------------------------------------------------
  /**
   * Query the component's ability to fire events
   * @param	void
   * @return	boolean	true if the component can fire events
   */
  public boolean isEnabled() {
     return _enabled;
  }
//-------------------------------------------------------------
  /**
   * Initialises the model and GUI.
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
   * hook up the adaptors 
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
// accessor methods for the bean properties etc
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
  public class textToModelAdaptor
    implements ActionListener {
    ZoomModel model;
    JTextField control;
    double VAL;
    String textstr = "";
    String msg = "the zoom control's text field expects a number\n such as 100 or 100.5";

  /**
   * Constructor
   */

    public textToModelAdaptor(JTextField cntrl, ZoomModel mdl) {
      model = mdl;
      control = cntrl;
    }

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
  public class modelToTextAdaptor
    implements ChangeListener {
    ZoomModel model;
    JTextField view;
    String valstr;

  /**
   * Constructor
   */
    public modelToTextAdaptor(ZoomModel mdl, JTextField vw) {
      view = vw;
      model = mdl;
    }

    public modelToTextAdaptor() {
    }

  /**
   * Adaptor, changes the <b>view</b> of the text field when the model changes.
   * @param	void
   * @return	boolean	true if the component can fire events
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
  // allows this bean to fire events
  public class modelToThisAdaptor
		implements ChangeListener {

  /**
   * Constructor
   */
    public modelToThisAdaptor() {
    }

  /**
   * Adaptor, fires an event from the bean when the model changes.
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
  // keep track of all the changeListeners to this model
  protected EventListenerList changeListeners =
     new EventListenerList();

  /**
   * Adds a ChangeListener for events fired from this bean.
   * @param	ChangeListener x
   * @return	void
   */
  public void addChangeListener(ChangeListener x) {
     changeListeners.add (ChangeListener.class, x);

     // bring it up to date with current state
     x.stateChanged(new ChangeEvent(this));
  }

  /**
   * Removes a ChangeListener for events fired from this bean.
   * @param	ChangeListener x
   * @return	void
   */
  public void removeChangeListener(ChangeListener x) {
     changeListeners.remove (ChangeListener.class, x);
  }

  /**
   * Fires an event if the zoom model changes.
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
