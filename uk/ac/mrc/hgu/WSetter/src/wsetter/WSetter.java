package wsetter;

import wsetter.*;

import javax.swing.*;
import java.awt.*;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

/**
 * View / Controller class for  WSetter.
 * @author Nick Burton
 * @see WSetterModel
 * @see SliderRangeModel
 */
public class WSetter extends WSetterGUI
                     implements WSetterConstants, Serializable {

   private boolean _enabled;
//-------------------------------------------------------------
   public WSetterModel mod1 = null;
   public SliderRangeModel sliderModel = null;

   protected sliderToModelAdaptor_A SToM_A = null;
   protected sliderToModelAdaptor_B SToM_B = null;
   protected textToModelAdaptor TToM = null;
   protected modelToTextAdaptor MToT = null;
   protected modelToSliderAdaptor MToS = null;
   protected limitToSliderAdaptor LToS = null;
   protected modelToThisAdaptor MToThis = null;
   // model properties
   private boolean slidingEvents = true;
   private int _type;
   private static int maxTextlen = 6;
//-------------------------------------------------------------
  /**
   * Creates a WSetter with double values and sliding events enabled.
   */
  public WSetter() {
     this(DOUBLE, true);
  }

  /**
   * Creates a WSetter with double values and sliding events as specified.
   */
  public WSetter(boolean sliding) {
     this(DOUBLE, sliding);
  }

  /**
   * Creates a WSetter with the specified value type and sliding events enabled.
   */
  public WSetter(int type) {
     this(type, true);
  }

  /**
   * Creates a WSetter with the specified value type and sliding events.
   */
  public WSetter(int type, boolean sliding) {
    slidingEvents = sliding;
    _type = type;
    _enabled = false;
    try {
      initWSetter();
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }
//-------------------------------------------------------------
  /**
   * Allow or inhibit the WSetter component's ability to fire events
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
   * @return	boolean	true if the WSetter component can fire events
   */
  public boolean isEnabled() {
     return _enabled;
  }
//-------------------------------------------------------------
  /**
   * Initialises the WSetter model and GUI.
   * @param	void
   * @return	void
   */
  private void initWSetter() throws Exception {

    // the model for the slider
    sliderModel = new SliderRangeModel();
    _slider.setModel(sliderModel);
    // the main WSetter model
    mod1 = new WSetterModel (_type);
    setModel();

  } // initWSetter()
//-------------------------------------------------------------
  /**
   * Hooks up the adaptors after first removing anty existing ones.
   * @param	void
   * @return	void
   */
  protected void setModel() {

     //-----------------------------------------
     // hook the adapters in place
     // (after removing any old ones)
     //-----------------------------------------

     // CONTROL ADAPTORS
     if (slidingEvents == true) {
	if(SToM_A != null) {
	   _slider.removeChangeListener(SToM_A);
        }
	SToM_A = new sliderToModelAdaptor_A(_slider, mod1);
	_slider.addChangeListener(SToM_A);
     } else {
	if(SToM_B != null) {
	   _slider.removeChangeListener(SToM_B);
        }
	SToM_B = new sliderToModelAdaptor_B(_slider, mod1);
	_slider.addChangeListener(SToM_B);
     }

     if(TToM != null) {
	textf.removeActionListener(TToM);
     }
     TToM = new textToModelAdaptor(textf, mod1);
     textf.addActionListener(TToM);
//...............................
     // VIEW ADAPTORS
     if(MToT != null) {
	mod1.removeChangeListener(MToT);
     }
     MToT = new modelToTextAdaptor(mod1, textf);
     mod1.addChangeListener(MToT);

     if(MToS != null) {
	mod1.removeChangeListener(MToS);
     }
     MToS = new modelToSliderAdaptor(mod1, _slider);
     mod1.addChangeListener(MToS);
//...............................
     // LIMIT ADAPTORS
     if(LToS != null) {
	mod1.removeLimitListener(LToS);
     }
     LToS = new limitToSliderAdaptor(mod1, _slider);
     mod1.addLimitListener(LToS);
//...............................
     // MODEL ADAPTORS
     if(MToThis != null) {
	mod1.removeChangeListener(MToThis);
     }
     //MToThis = new modelToThisAdaptor(mod1, this);
     MToThis = new modelToThisAdaptor(this);
     mod1.addChangeListener(MToThis);

  } // setModel

  public int getType() {
     return mod1.getType();
  }

//===========================================
// accessor methods for the bean properties etc
//===========================================
   public Vector getValue() {
      return mod1.getValue();
   }
   public Vector getMin() {
      return mod1.getMin();
   }
   public Vector getMax() {
      return mod1.getMax();
   }
//.....................................
   public int getWidth() {
      return _width;
   }
   public int getHeight() {
      return _height;
   }
   public int getLabelWidth() {
      return labelWidth;
   }
   public int getTextWidth() {
      return textWidth;
   }
   public Color getBgc() {
      return bgc;
   }
   public boolean getSlidingEvents() {
      return slidingEvents;
   }
   public String getSliderLabel() {
      return sliderLabel;
   }
   public boolean isAdjusting() {
      return _slider.getValueIsAdjusting();
   }
//.....................................
   public void setValue(double val) {
      mod1.setValue(val);
   }
   public void setMin(double val) {
      mod1.setMin(val);
   }
   public void setMax(double val) {
      mod1.setMax(val);
   }
//.....................................
   public void setWidth(int w) {
      _width = w;
      if(_width < _minWidth) _width = _minWidth;
      initGUI();
   }
   public void setHeight(int h) {
      _height = h;
      if(_height < _minHeight) _height = _minHeight;
      initGUI();
   }
   public void setLabelWidth(int w) {
      labelWidth = w;
      initGUI();
   }
   public void setTextWidth(int w) {
      textWidth = w;
      initGUI();
   }
   public void setBgc(Color col) {
      bgc = col;
      initGUI();
   }
   /*
   public  void setSlidingEvents(boolean bool) {
      slidingEvents = bool;
      setModel();
   }
   */
   public void setSliderLabel(String labl) {
      sliderLabel = labl;
      initGUI();
   }
//.....................................
   public void setSliderEnabled(boolean bool) {
      _slider.setEnabled(bool);
   }

   public boolean isSliderEnabled() {
      return _slider.isEnabled();
   }

//===========================================
// inner classes
//===========================================
// adapters for MVC model

// CONTROL ADAPTORS
//---------------------------------------
  /**
   * Handles continuous changes in slider position
   */
  public static class sliderToModelAdaptor_A
    implements ChangeListener {

    WSetterModel model;
    JSlider control;
    SliderRangeModel sliderMod;

  /**
   */
    public sliderToModelAdaptor_A(JSlider cntrl, WSetterModel mdl) {
      model = mdl;
      control = cntrl;
      sliderMod = (SliderRangeModel) control.getModel();
    }

  /**
   * Event handler, for use when slider is moving.
   * @param	ChangeEvent e
   * @return	void
   */
    public void stateChanged(ChangeEvent e) {
      model.setValue(sliderMod.getDoubleValue());
    }
  } // class sliderToModelAdaptor

//---------------------------------------
  /**
   * Handles a step change in slider position
   *
   */
  public static class sliderToModelAdaptor_B
    implements ChangeListener {
    WSetterModel model;
    JSlider control;
    SliderRangeModel sliderMod;

  /**
   */
    public sliderToModelAdaptor_B(JSlider cntrl, WSetterModel mdl) {
      model = mdl;
      control = cntrl;
      sliderMod = (SliderRangeModel) control.getModel();
    }

  /**
   * Event handler, for use when slider has stopped. 
   * @param	ChangeEvent e
   * @return	void
   */
    public void stateChanged(ChangeEvent e) {
      // only do it when the slider isn't moving
      if (!control.getValueIsAdjusting()) {
	 model.setValue(sliderMod.getDoubleValue());
      }
    }
  } // class sliderToModelAdaptor

//---------------------------------------
  /**
   * Updates the WSetter model when a numeric value is entered into the WSetter's text field.
   *
   */
  public static class textToModelAdaptor
    implements ActionListener {
    WSetterModel model;
    JTextField control;
    double VAL;
    String textstr = "";
    String msg = "the slider's text field expects a number\n such as 10 or 10.5";

  /**
   */
    public textToModelAdaptor(JTextField cntrl, WSetterModel mdl) {
      model = mdl;
      control = cntrl;
    }

  /**
   * Event handler, sets the model to the value in the text field.
   * @param	ActionEvent e
   * @return	void
   */
    public void actionPerformed(ActionEvent e) {
       try {
          textstr = control.getText();
	  VAL = Double.parseDouble(textstr);
	  model.setValue(VAL);
       }
       catch(NullPointerException npe) {
          System.out.println(
	       "null pointer exception in WSetter textToModelAdaptor");
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
   * Updates the slider view of the WSetter when it's model has changed.
   *
   */
  public static class modelToSliderAdaptor
    implements ChangeListener {
    WSetterModel model;
    JSlider view;
    SliderRangeModel sliderMod;
    Vector vec = null;
    int type;

  /**
   */
    public modelToSliderAdaptor(WSetterModel mdl, JSlider vw) {
      view = vw;
      model = mdl;
      sliderMod = (SliderRangeModel)view.getModel();
      type = model.getType();
    }

  /**
   * Adaptor, updates the slider view of the WSetter when it's model has changed.
   * @param	ChangeEvent e
   * @return	void
   */
    public void stateChanged(ChangeEvent e) {
      vec = model.getValue();
      switch(type) {
         case INTEGER:
	    Integer iVal = (Integer)vec.elementAt(0);
	    sliderMod.setDoubleValue((double)iVal.intValue());
	    break;
         case FLOAT:
	    Float fVal = (Float)vec.elementAt(0);
	    sliderMod.setDoubleValue((double)fVal.floatValue());
	    break;
         case DOUBLE:
	    Double dVal = (Double)vec.elementAt(0);
	    sliderMod.setDoubleValue(dVal.doubleValue());
	    break;
	 default:
	    break;
      }
      view.revalidate();
    }
  } // class modelToSliderAdaptor

//---------------------------------------
  /**
   * Updates the text field view of the WSetter when it's model has changed.
   *
   */
  public static class modelToTextAdaptor
    implements ChangeListener {
    WSetterModel model;
    JTextField view;
    int type;
    Vector vec = null;

  /**
   */
    public modelToTextAdaptor(WSetterModel mdl, JTextField vw) {
      view = vw;
      model = mdl;
      type = model.getType();
    }

    String valstr;
  /**
   * Adaptor, updates the text field view of the WSetter when it's model has changed.
   * @param	ChangeEvent e
   * @return	void
   */
    public void stateChanged(ChangeEvent e) {
      vec = model.getValue();
      switch(type) {
         case INTEGER:
	    Integer iVal = (Integer)vec.elementAt(0);
	    valstr = (iVal.toString());
	    break;
         case FLOAT:
	    Float fVal = (Float)vec.elementAt(0);
	    valstr = (fVal.toString());
	    break;
         case DOUBLE:
	    Double dVal = (Double)vec.elementAt(0);
	    valstr = (dVal.toString());
	    break;
	 default:
	    break;
      }

      if(valstr.length() > maxTextlen) {
	 view.setText(valstr.substring(0,maxTextlen));
      } else {
	 view.setText(valstr);
      }
    }
  } // class modelToTextAdaptor

//---------------------------------------
// LIMIT ADAPTORS
// changes slider when WSetter max/min limits are changed
//---------------------------------------
  /**
   * Updates the WSetter's slider's max/min limits when the WSetter model's limits have changed.
   *
   */
  public static class limitToSliderAdaptor
    implements LimitListener {
    WSetterModel model;
    JSlider view;
    SliderRangeModel sliderMod;
    Vector maxVec = null;
    Vector minVec = null;
    int type;

  /**
   */
    public limitToSliderAdaptor(WSetterModel mdl, JSlider vw) {
      view = vw;
      model = mdl;
      sliderMod = (SliderRangeModel)view.getModel();
      type = model.getType();
    }

  /**
   * Adaptor, updates the slider's max/min limits when the WSetter model's limits have changed.
   */
    public void limitChanged(LimitEvent e) {
      maxVec = model.getMax();
      minVec = model.getMin();
      switch(type) {
         case INTEGER:
	    Integer imaxVal = (Integer)maxVec.elementAt(0);
	    Integer iminVal = (Integer)minVec.elementAt(0);
	    sliderMod.setMaximum(imaxVal.intValue());
	    sliderMod.setMinimum(iminVal.intValue());
	    break;
         case FLOAT:
	    Float fmaxVal = (Float)maxVec.elementAt(0);
	    Float fminVal = (Float)minVec.elementAt(0);
	    sliderMod.setMaximum((int)(fmaxVal.floatValue()));
	    sliderMod.setMinimum((int)(fminVal.floatValue()));
	    break;
         case DOUBLE:
	    Double dmaxVal = (Double)maxVec.elementAt(0);
	    Double dminVal = (Double)minVec.elementAt(0);
	    sliderMod.setMaximum((int)(dmaxVal.intValue()));
	    sliderMod.setMinimum((int)(dminVal.intValue()));
	    break;
	 default:
	    break;
      }
      view.revalidate();
    }
  } // class limitToSliderAdaptor

//---------------------------------------
  // MODEL ADAPTORS
//---------------------------------------
  /**
   * Fires an event when the WSetter's model has changed.
   *
   */
  // allows this bean to fire events
  public static class modelToThisAdaptor
		implements ChangeListener {

    WSetter target;
  /**
   */
    public modelToThisAdaptor(WSetter wst) {
      target = wst;
    }

  /**
   * Adaptor, fires an event when the WSetter's model has changed.
   * @param	ChangeEvent e
   * @return	void
   */
    public void stateChanged(ChangeEvent e) {
      target.fireChange();
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
   * Fires an event if the slider model changes.
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

} // class WSetter
