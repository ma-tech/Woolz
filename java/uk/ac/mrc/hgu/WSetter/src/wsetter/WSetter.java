package wsetter;

import wsetter.*;

import javax.swing.*;
import java.awt.*;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;
//import com.borland.jbcl.layout.*;

public class WSetter extends JPanel
                     implements Serializable {


  public WSetter() {
    try {
      jbInit();
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }

// adapters for MVC model

// CONTROL ADAPTORS
//---------------------------------------
  public static class sliderToModelAdaptor_A
    implements ChangeListener {
    WlzFltModel model;
    JSlider control;
    SliderRangeModel sliderMod;

    public sliderToModelAdaptor_A(JSlider cntrl, WlzFltModel mdl) {
      model = mdl;
      control = cntrl;
    }

    public void stateChanged(ChangeEvent e) {
      sliderMod = (SliderRangeModel) control.getModel();
      model.setValue((float)(sliderMod.getDoubleValue()));
    }
  } // class sliderToModelAdaptor

//---------------------------------------
  public static class sliderToModelAdaptor_B
    implements ChangeListener {
    WlzFltModel model;
    JSlider control;
    SliderRangeModel sliderMod;

    public sliderToModelAdaptor_B(JSlider cntrl, WlzFltModel mdl) {
      model = mdl;
      control = cntrl;
    }

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
    public textToModelAdaptor(JTextField cntrl, WlzFltModel mdl) {
      model = mdl;
      control = cntrl;
    }

    float val = 0.0f;
    Float VAL;
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

    public modelToSliderAdaptor(WlzFltModel mdl, JSlider vw) {
      view = vw;
      model = mdl;
    }

    public void stateChanged(ChangeEvent e) {
      // get the slider's model
      sliderMod = (SliderRangeModel)view.getModel();
      sliderMod.setDoubleValue((double)model.getValue());
      view.invalidate();
    }
  } // class modelToSliderAdaptor

//---------------------------------------
  public static class modelToTextAdaptor
    implements ChangeListener {
    WlzFltModel model;
    JTextField view;
    public modelToTextAdaptor(WlzFltModel mdl, JTextField vw) {
      view = vw;
      model = mdl;
    }

    Float flote;
    public void stateChanged(ChangeEvent e) {
      flote = new Float(model.getValue());
      view.setText(flote.toString());
      view.invalidate();
    }
  } // class modelToTextAdaptor

//---------------------------------------
  // MODEL ADAPTORS
//---------------------------------------
  // allows this bean to fire events
  public static class modelToThisAdaptor
    implements ChangeListener {
    WlzFltModel source;
    WSetter target;
    public modelToThisAdaptor(WlzFltModel mdl1, WSetter wst) {
      source = mdl1;
      target = wst;
    }

    public void stateChanged(ChangeEvent e) {
      target.fireChange();
    }
  } // class modelToThisAdaptor


//*******************************************
  // GUI stuff
  JPanel textPanel = new JPanel();
  JTextField textf = new JTextField();

  JPanel sliderPanel = new JPanel();
  JSlider _slider = new JSlider();

  JPanel labelSliderTextPanel = new JPanel();

  JPanel labelPanel = new JPanel();
  JLabel paramLabel = new JLabel("", SwingConstants.CENTER);
//*******************************************
//===========================================
  private void jbInit() throws Exception {

    sliderModel = new SliderRangeModel();
    _slider.setModel(sliderModel);
    setGUI();
    setModel();

  } // jbInit()
//===========================================

  protected void setGUI() {
    sliderH = height;
    sliderW = width - (labelWidth+textWidth);

    this.setLayout(new BorderLayout());
    //this.setLayout(new FlowLayout());
    this.setBackground(internalBgc);
    //this.setPreferredSize(new Dimension(width+pad, height+pad));
    this.setPreferredSize(new Dimension(0, height+pad));
    //this.setMaximumSize(new Dimension(width+pad, height+pad));

    textf.setBackground(bgc);
    textf.setPreferredSize(new Dimension(textWidth, sliderH+pad));

    textPanel.setBackground(internalBgc);
    textPanel.setLayout(new BorderLayout(0, pad));
    textPanel.setPreferredSize(new Dimension(textWidth+pad, sliderH+pad));
    textPanel.add(textf, BorderLayout.NORTH);

    //slider.setLayout(new BorderLayout());
    _slider.setLayout(new FlowLayout());
    _slider.setPreferredSize(new Dimension(sliderW, sliderH));
    //slider.setExtent(sliderExtent);
    sliderPanel.setBackground(internalBgc);
    sliderPanel.setLayout(new BorderLayout());
    //sliderPanel.setPreferredSize(new Dimension(sliderW, sliderH));
    sliderPanel.setPreferredSize(new Dimension(0, sliderH));
    sliderPanel.add(_slider, BorderLayout.CENTER);

    paramLabel.setBackground(internalBgc);
    paramLabel.setPreferredSize(new Dimension(labelWidth, sliderH));
    paramLabel.setText(sliderLabel);

    labelPanel.setBackground(bgc);
    labelPanel.setLayout(new BorderLayout(pad, pad));
    labelPanel.setPreferredSize(new Dimension(labelWidth, sliderH));
    labelPanel.add(paramLabel, BorderLayout.NORTH);

    labelSliderTextPanel.setBackground(internalBgc);
    labelSliderTextPanel.setLayout(new BorderLayout());
    labelSliderTextPanel.setPreferredSize(
			 new Dimension(width, height));
    labelSliderTextPanel.add(labelPanel, BorderLayout.WEST);
    labelSliderTextPanel.add(sliderPanel, BorderLayout.CENTER);
    labelSliderTextPanel.add(textPanel, BorderLayout.EAST);
    this.add(labelSliderTextPanel, BorderLayout.NORTH);

  } // setGUI

//-------------------------------------------------------------
  // if the programmer wants to change model max & min
  // we need to set up the adaptors again
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
     MToThis = new modelToThisAdaptor(mod1, this);
     mod1.addChangeListener(MToThis);

  } // setModel

//-------------------------------------------------------------
  // we need to generate change events
  // whenever the WlzFltModel changes

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
  } // fireChange

//-------------------------------------------------------------
// accessor methods for the bean properties etc
   public float getValue() {
      return mod1.getValue();
   }
//.....................................
   public int getWidth() {
      return width;
   }
   public int getHeight() {
      return height;
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
      width = w;
      if(width < minWidth) width = minWidth;
      setGUI();
   }
   public void setHeight(int h) {
      height = h;
      if(height < minHeight) height = minHeight;
      setGUI();
   }
   public void setLabelWidth(int w) {
      labelWidth = w;
      setGUI();
   }
   public void setTextWidth(int w) {
      textWidth = w;
      setGUI();
   }
   public void setBgc(Color col) {
      bgc = col;
      setGUI();
   }
   public  void setFloatingPoint(boolean bool) {
      floatingPoint = bool;
   }
   public  void setSlidingEvents(boolean bool) {
      slidingEvents = bool;
   }
   public void setSliderLabel(String labl) {
      sliderLabel = labl;
      setGUI();
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
   public void setEnabled(boolean bool) {
      _slider.setEnabled(bool);
   }

   //public void setModelType(WlzFltMdl.WlzType WlzFltMdl.WlzType.FLOAT)

//-------------------------------------------------------------
   public WlzFltModel mod1;
   public SliderRangeModel sliderModel;
   //....................................
   protected sliderToModelAdaptor_A SToM_A = null;
   protected sliderToModelAdaptor_B SToM_B = null;
   protected textToModelAdaptor TToM = null;
   protected modelToTextAdaptor MToT = null;
   protected modelToSliderAdaptor MToS = null;
   protected modelToThisAdaptor MToThis = null;

//-------------------------------------------------------------
   // GUI properties
   // private int width = 250;
   private int width = 800;
   private int height = 15;
   private int textWidth = 50;
   private int labelWidth = 40;
   private int pad = 1;
   private int minWidth = 150; // not accessible
   private int minHeight = 10; // not accessible

   private Color bgc = new Color(230, 230, 230);
   private Color internalBgc = new Color(200, 200, 200); // not accessible

   private boolean floatingPoint = true; // if false => int
   private boolean slidingEvents = true;

   private String sliderLabel = "";
   //....................................
   // model properties
   private WlzFltModel.WlzType modelType = WlzFltModel.FLOAT;
   private float modelMin = 0.0f;
   private float modelMax = 1000.0f;
   private float modelInitVal = 0.0f;
   //....................................
   // constants to define GUI shape
   int sliderH;
   int sliderW;

} // class WSetter
