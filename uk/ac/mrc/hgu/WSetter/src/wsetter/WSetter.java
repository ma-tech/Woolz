package wsetter;

import wsetter.*;

import javax.swing.*;
import java.awt.*;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;

/**
 * View / Controller class for slider bean.
 * <br>Uses the <b>Model View Controller</b> paradigm.
 * @author Nick Burton
 * @see WlzFltModel
 * @see SliderRangeModel
 */
public class WSetter extends WSetterGUI
                     implements Serializable {

   private boolean _enabled;
//-------------------------------------------------------------
   public WlzFltModel mod1;
   public SliderRangeModel sliderModel;

   protected sliderToModelAdaptor_A SToM_A = null;
   protected sliderToModelAdaptor_B SToM_B = null;
   protected textToModelAdaptor TToM = null;
   protected modelToTextAdaptor MToT = null;
   protected modelToSliderAdaptor MToS = null;
   protected modelToThisAdaptor MToThis = null;
   // model properties
   private boolean floatingPoint = true; // if false => int
   private boolean slidingEvents = true;
   private WlzFltModel.WlzType modelType = WlzFltModel.FLOAT;
   private float modelMin = 0.0f;
   private float modelMax = 1000.0f;
   private float modelInitVal = 0.0f;
//-------------------------------------------------------------
  /**
   * Constructor
   */
  public WSetter() {
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
   * Initialises the slider model and GUI.
   * @param	void
   * @return	void
   */
  private void initWSetter() throws Exception {

    sliderModel = new SliderRangeModel();
    _slider.setModel(sliderModel);
    setModel();

  } // initWSetter()
//-------------------------------------------------------------
  // if the programmer wants to change model max & min
  // we need to set up the adaptors again
  /**
   * Resets the model to new max and min values.
   * @param	void
   * @return	void
   */
  protected void setModel() {

     // the main model
     mod1 = null;
     mod1 = new WlzFltModel (modelType,
			     modelMin,
			     modelMax);

     mod1.setValue(modelInitVal);

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

     // MODEL ADAPTORS
     if(MToThis != null) {
	mod1.removeChangeListener(MToThis);
     }
     //MToThis = new modelToThisAdaptor(mod1, this);
     MToThis = new modelToThisAdaptor(this);
     mod1.addChangeListener(MToThis);

  } // setModel

//===========================================
// accessor methods for the bean properties etc
//===========================================
   public float getValue() {
      return mod1.getValue();
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
   public boolean getFloatingPoint() {
      return floatingPoint;
   }
   public boolean getSlidingEvents() {
      return slidingEvents;
   }
   public String getSliderLabel() {
      return sliderLabel;
   }
//.....................................
   public float getModelMin() {
      return modelMin;
   }
   public float getModelMax() {
      return modelMax;
   }
   public float getModelInitVal() {
      return modelInitVal;
   }

//.....................................
   public void setValue(float val) {
      mod1.setValue(val);
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
   public  void setFloatingPoint(boolean bool) {
      floatingPoint = bool;
   }
   public  void setSlidingEvents(boolean bool) {
      slidingEvents = bool;
   }
   public void setSliderLabel(String labl) {
      sliderLabel = labl;
      initGUI();
   }
//.....................................
   public void setModelMin(float mmin) {
      modelMin = mmin;
      setModel();
      sliderModel.setMinimum((int)mod1.getMinimum());
   }
   public void setModelMax(float mmax) {
      modelMax = mmax;
      setModel();
      sliderModel.setMaximum((int)mod1.getMaximum());
   }
   public void setModelInitVal(float minit) {
      modelInitVal = minit;
      setModel();
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
   * Inner class to handle continuous changes in slider position
   */
  public static class sliderToModelAdaptor_A
    implements ChangeListener {

    WlzFltModel model;
    JSlider control;
    SliderRangeModel sliderMod;

  /**
   * Constructor
   */
    public sliderToModelAdaptor_A(JSlider cntrl, WlzFltModel mdl) {
      model = mdl;
      control = cntrl;
    }

  /**
   * Event handler, for use when slider has stopped.
   * @param	ChangeEvent e
   * @return	void
   */
    public void stateChanged(ChangeEvent e) {
      sliderMod = (SliderRangeModel) control.getModel();
      model.setValue((float)(sliderMod.getDoubleValue()));
    }
  } // class sliderToModelAdaptor

//---------------------------------------
  /**
   * Inner class to handle a step change in slider position
   *
   */
  public static class sliderToModelAdaptor_B
    implements ChangeListener {
    WlzFltModel model;
    JSlider control;
    SliderRangeModel sliderMod;

  /**
   * Constructor
   */
    public sliderToModelAdaptor_B(JSlider cntrl, WlzFltModel mdl) {
      model = mdl;
      control = cntrl;
    }

  /**
   * Event handler, for use when slider has is moving. 
   * @param	void
   * @return	boolean	true if the component can fire events
   */
    public void stateChanged(ChangeEvent e) {
      // only do it when the slider isn't moving
      if (!control.getValueIsAdjusting()) {
	 // get the model
	 sliderMod = (SliderRangeModel) control.getModel();
	 model.setValue((float)(sliderMod.getDoubleValue()));
      }
    }
  } // class sliderToModelAdaptor

//---------------------------------------
  public static class textToModelAdaptor
    implements ActionListener {
    WlzFltModel model;
    JTextField control;
  /**
   * Constructor
   */
    public textToModelAdaptor(JTextField cntrl, WlzFltModel mdl) {
      model = mdl;
      control = cntrl;
    }

    float val = 0.0f;
    Float VAL;
  /**
   * Event handler, sets the slider to the value in the text field.
   * @param	ActionEvent e
   * @return	void
   */
    public void actionPerformed(ActionEvent e) {
       VAL = new Float(control.getText());
       val = VAL.floatValue();
       model.setValue(val);
    }
  } // class textToModelAdaptor

//---------------------------------------
// VIEW ADAPTORS
//---------------------------------------
  public static class modelToSliderAdaptor
    implements ChangeListener {
    WlzFltModel model;
    JSlider view;
    SliderRangeModel sliderMod;

  /**
   * Constructor
   */
    public modelToSliderAdaptor(WlzFltModel mdl, JSlider vw) {
      view = vw;
      model = mdl;
    }

  /**
   * Adaptor, changes the <b>view</b> of the slider when the model changes.
   * @param	void
   * @return	boolean	true if the component can fire events
   */
    public void stateChanged(ChangeEvent e) {
      // get the slider's model
      sliderMod = (SliderRangeModel)view.getModel();
      sliderMod.setDoubleValue((double)model.getValue());
      view.revalidate();
    }
  } // class modelToSliderAdaptor

//---------------------------------------
  public static class modelToTextAdaptor
    implements ChangeListener {
    WlzFltModel model;
    JTextField view;
  /**
   * Constructor
   */
    public modelToTextAdaptor(WlzFltModel mdl, JTextField vw) {
      view = vw;
      model = mdl;
    }

    Float flote;
  /**
   * Adaptor, changes the <b>view</b> of the text field when the model changes.
   * @param	void
   * @return	boolean	true if the component can fire events
   */
    public void stateChanged(ChangeEvent e) {
      flote = new Float(model.getValue());
      view.setText(flote.toString());
    }
  } // class modelToTextAdaptor

//---------------------------------------
  // MODEL ADAPTORS
//---------------------------------------
  // allows this bean to fire events
  public static class modelToThisAdaptor
		implements ChangeListener {

    WSetter target;
  /**
   * Constructor
   */
    public modelToThisAdaptor(WSetter wst) {
      target = wst;
    }

  /**
   * Adaptor, fires an event from the bean when the model changes.
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
  // keep track of all the listeners to this model
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
