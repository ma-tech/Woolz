package sectionViewer;

import sectionViewer.*;

import java.lang.reflect.Method;
import java.lang.reflect.InvocationTargetException;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;
import javax.help.*;
import java.net.*;
import java.util.*;

import com.sun.image.codec.jpeg.*;

import wsetter.*;
import zoom.*;
import uk.ac.mrc.hgu.Wlz.*;

/** 
 *   Combined <b>View</b> & <b>Controller</b> class for SectionViewer component.
 */
public class SectionViewer
    extends SectionViewerGUI
    implements Serializable, WlzObjectType, WSetterConstants {

  /**  Toggle for display of debugging information. */
  private final boolean _debug = false;

  /**  True if SectionViewer can fire events */
  private boolean _enabled = true;

  //protected boolean _parentIsRoot = true;

  /**
   *   Instance of class that implements SVParent interface.
   *   <br>SVParent must be implemented by the class responsible for 
   *   creating SectionViewers.
   */
  public Object _parent = null;

  /**   Cursor used for operations within a grey-level image. */
  protected Cursor xhairCursor = new Cursor(Cursor.CROSSHAIR_CURSOR);

  /**   The default cursor shape. */
  protected Cursor defCursor = Cursor.getDefaultCursor();

  /**   Model that encapsulates Woolz object data. */
  protected WlzObjModel _OBJModel = null; // grey level

  /**   Model that encapsulates View Structure data. */
  protected ViewStructModel _VSModel = null;

  /**   Instance of class responsible for drawing the screen image. */
  public WlzImgView _imgV = null;

  /**
   *   Prevents multiple events from fixed point
   *   <em>PointEntry</em> dialogue.
   */
  private int _FPBounce = 0; // stops initial false event from PointEntry

  /**
   *   Prevents multiple events from 2nd fixed point
   *   <em>PointEntry</em> dialogue.
   */
  private int _APBounce = 0;

  /*  Unused variable.
  private WlzObject _anatomyObj = null; */
  
  /**   Woolz object representing a threshold constraint region */
  private WlzObject _constraintObj = null;

  /**   Instance of class responsible for building anatomy menus. */
  private AnatomyBuilder _anatBuilder = null;

  /**  
   *   Collection of full pathnames to anatomy components
   *   in the <em>embryo</em> hierarchy.
   */
  private Stack _embFileStack = null;

  /**  
   *   Collection of <em>Woolz</em> objects representing anatomy components
   *   in the <em>embryo</em> hierarchy.
   */
  private Stack _embObjStack = null;

  /**  
   *   Collection of full pathnames to anatomy components
   *   in the <em>extra-embryonic component</em> hierarchy.
   */
  private Stack _xembFileStack = null;

  /**  
   *   Collection of <em>Woolz</em> objects representing anatomy components
   *   in the <em>extra-embryonic component</em> hierarchy.
   */
  private Stack _xembObjStack = null;

  /**  
   *   Collection of <em>File</em> objects representing
   *   potential <em>mouse-click anatomy</em> components.
   */
  private Stack _maybeFil = null;

  /**  
   *   Collection of <em>Woolz</em> objects representing
   *   potential <em>mouse-click anatomy</em> components.
   */
  private Stack _maybeObj = null;

  /**
   *   Description of initial SectionViewer orientation.
   *   <br>Possible descriptions are:
   *   <ul>
   *      <li>XY,  (corresponds to Pitch = 0, Yaw = 0, Roll = 0)</li>
   *      <li>YZ,  (corresponds to Pitch = 90, Yaw = 0, Roll = 90)</li>
   *      <li>ZX,  (corresponds to Pitch = 90, Yaw = 90, Roll = 90)</li>
   *   </ul>
   */
  protected String viewType;

  /**   Full pathname of <em>mouse-click anatomy</em> component. */
  protected String _anatomyStr;

  /**   True if <em>mouse-click anatomy</em> is enabled. */
  private boolean _anatomy; // for mouse click anatomy

  /**   True if <em>thresholding</em> is enabled. */
  private boolean _thresholding;

  /**   True if <em>fixed line rotation</em> is enabled. */
  private boolean _fixedLineRotation;

  /**   True if <em>threshold constraint</em> is enabled. */
  private boolean _threshConstraint;

  /** 
   *   True if <em>fixed point</em> has been selected
   *   from the <em>show</em> menu.
   */
  protected boolean _fixedPoint = false;

  /**
   *   True if <em>fixed line</em> has been selected
   *   from the <em>show</em> menu.
   */
  private boolean _axisPoint = false;

  /**
   *   True if <em>change fixed point using mouse button</em>
   *   has been selected from the
   *   <em>control / fixed point</em> menu.
   */
  private boolean _setFixedPoint = false;

  /**
   *   True if <em>change fixed line using mouse button</em>
   *   has been selected from the
   *   <em>control / fixed line end-point</em> menu.
   */
  private boolean _setAxisPoint = false;

  /**
   *   True if <em>yaw pitch roll</em>
   *   has been selected from the
   *   <em>control / rotation</em> menu.
   */
  private boolean _showStdCntrls = false;

  /**
   *   True if <em>fixed line</em>
   *   has been selected from the
   *   <em>control / rotation</em> menu.
   */
  private boolean _showUsrCntrls = false;

  /**
   *   True if <em>intersection lines</em> has been selected
   *   from the <em>show</em> menu.
   */
  private boolean _showIntersection = false;

  /**
   *   True if <em>openView()</em> has been called
   *   but has not completed.
   */
  private boolean _openingView = false;

  /**
   *   Flag toggled when the <em>invert button</em> is pressed.
   *   <br>If true the grey-level image has its grey values inverted.
   */
  private boolean _inverted = false;

  /**   flag required by <em>tie point</em> application. */
  public boolean _drawfixLine = false;

  /**
   *   Event handler for <em>file</em> menu selections.
   */
  fileMenuHandler handler_1 = null;

  /**
   *   Event handler for <em>control</em> menu selections.
   */
  controlMenuHandler handler_2 = null;

  /**
   *   Event handler for <em>show</em> menu selections.
   */
  showMenuHandler handler_3 = null;

  /**
   *   Event handler for <em>help</em> menu selections.
   */
  helpMenuHandler handler_4 = null;

  // temporarily disabled ... don't remove
  /**
   *   Event handler for <em>threshold</em> menu selections.
   *   <p>Not in use.
   */
  thresholdMenuHandler handler_5 = null;

  //--------------------------------------

  /**
   *   Event handler for <em>secColorClt</em>
   *   (colourChooser button) control.
   */
  public planeColChooser colourHandler = null;

  /**
   *   Event handler for <em>invertButton</em> control.
   */
  public invertButtonHandler invertHandler = null;

  /**   System file separator ("/" or "\"). */
  private String SLASH = System.getProperty("file.separator");

  /**
   *   Container for SectionViewer component.
   *   <br>JFrame for external components.
   *   <br>JInternalFrame for internal components.
   */
  protected Container _frame = null;

  //=========================================================
  // constructor
  //=========================================================

  /**
   *   Constructs a SectionViewer with the given
   *   initial orientation and parent.
   *   @param viewstr the initial orientation.
   *   @param parent instance of the class that created this SectionViewer.
   *   <br>parent (<b>Must</b> implement <em>SVParent</em> interface.)
   */
  public SectionViewer(String viewstr, Object parent) {

    if (_debug) System.out.println("enter SectionViewer");

    _parent = parent;

    _anatomy = false;
    _thresholding = false;
    _threshConstraint = false;
    _fixedLineRotation = false;

    setCursor(defCursor);
    addFeedback();

    // instantiate the event handlers
    handler_1 = new fileMenuHandler();
    handler_2 = new controlMenuHandler();
    handler_3 = new showMenuHandler();
    handler_4 = new helpMenuHandler();
    handler_5 = new thresholdMenuHandler();

    // attach event handlers to menu items
    fileMenu_1.addActionListener(handler_1);
    fileMenu_2.addActionListener(handler_1);
    fileMenu_3.addActionListener(handler_1);
    fileMenu_4.addActionListener(handler_1);

    controlMenu_1_1.addActionListener(handler_2);
    controlMenu_1_2.addActionListener(handler_2);
    controlMenu_1_3.addActionListener(handler_2);
    controlMenu_2_1.addActionListener(handler_2);
    controlMenu_2_2.addActionListener(handler_2);
    controlMenu_2_3.addActionListener(handler_2);
    controlMenu_3_1.addActionListener(handler_2);
    controlMenu_3_2.addActionListener(handler_2);
    controlMenu_3_3.addActionListener(handler_2);
    controlMenu_4_1.addActionListener(handler_2);
    controlMenu_4_2.addActionListener(handler_2);
    controlMenu_4_3.addActionListener(handler_2);
    controlMenu_5.addActionListener(handler_2);

    showMenu_1.addActionListener(handler_3);
    showMenu_2.addActionListener(handler_3);
    showMenu_3.addActionListener(handler_3);
    showMenu_4.addActionListener(handler_3);
    showMenu_5.addActionListener(handler_3);

    thresholdMenu_1.addActionListener(handler_5);
    thresholdMenu_2.addActionListener(handler_5);
    thresholdMenu_3.addActionListener(handler_5);
    thresholdMenu_4.addActionListener(handler_5);

    helpMenu_1.addActionListener(handler_4);
    helpMenu_2.addActionListener(handler_4);
    helpMenu_3.addActionListener(handler_4);

    colourHandler = new planeColChooser();
    secColorClt.addActionListener(colourHandler);

    invertHandler = new invertButtonHandler();
    invertButton.addActionListener(invertHandler);

    zoomSetter.setEnabled(true);
    distSetter.setEnabled(true);
    pitchSetter.setEnabled(true);
    yawSetter.setEnabled(true);
    rollSetter.setEnabled(true);
    rotSetter.setEnabled(true);

    viewType = viewstr.substring(0, 2);
    if (viewType.equals("XY")) {
      if (pitchSetter != null) {
        pitchSetter.setValue(0.0);
      }
      if (yawSetter != null) {
        yawSetter.setValue(0.0);
      }
      if (rollSetter != null) {
        rollSetter.setValue(0.0);
      }
    }
    else if (viewType.equals("YZ")) {
      if (pitchSetter != null) {
        pitchSetter.setValue(90.0);
      }
      if (yawSetter != null) {
        yawSetter.setValue(0.0);
      }
      if (rollSetter != null) {
        rollSetter.setValue(90.0);
      }
    }
    else if (viewType.equals("ZX")) {
      if (pitchSetter != null) {
        pitchSetter.setValue(90.0);
      }
      if (yawSetter != null) {
        yawSetter.setValue(90.0);
      }
      if (rollSetter != null) {
        rollSetter.setValue(90.0);
      }
    }
    rollSetter.setSliderEnabled(false);
    if (rotSetter != null) {
      rotSetter.setValue(0.0);
      rotSetter.setSliderEnabled(false);
    }
    resetFeedbackText();

    if (_debug)
      System.out.println("exit SectionViewer");
  } // constructor

//-------------------------------------------------------------
// methods for adding / removing panels from the gui
//-------------------------------------------------------------
  /**
   *   Makes <em>yaw pitch</em> and <em>roll</em> controls visible.
   */
  protected void addStandardControls() {
    Dimension dim;
    if (_showUsrCntrls) {
      dim = new Dimension(totalW, transH);
    }
    else {
      dim = new Dimension(totalW, pyrH);
    }
    transientPanel.setPreferredSize(dim);
    transientPanel.setMinimumSize(dim);
    transientPanel.add(pitchYawRollPanel, BorderLayout.NORTH);
    revalidate();
    _bigPanel.repaint();
  }

  /**
   *   Makes <em>yaw pitch</em> and <em>roll</em> controls invisible.
   */
  protected void removeStandardControls() {
    Dimension dim;
    if (_showUsrCntrls) {
      dim = new Dimension(totalW, rotH);
    }
    else {
      dim = new Dimension(totalW, pad);
    }
    transientPanel.setPreferredSize(dim);
    transientPanel.setMinimumSize(dim);
    transientPanel.remove(pitchYawRollPanel);
    transientPanel.repaint();
    revalidate();
    _bigPanel.repaint();
  }

  //...............................

  /**
   *   Makes <em>fixed line</em> control visible.
   */
  protected void addUserControls() {
    Dimension dim;
    if (_showStdCntrls) {
      dim = new Dimension(totalW, transH);
    }
    else {
      dim = new Dimension(totalW, rotH);
    }
    transientPanel.setPreferredSize(dim);
    transientPanel.setMinimumSize(dim);
    transientPanel.add(userDefinedRotPanel, BorderLayout.SOUTH);
    revalidate();
    _bigPanel.repaint();

  }

  /**
   *   Makes <em>fixed line</em> control invisible.
   */
  protected void removeUserControls() {
    Dimension dim;
    if (_showStdCntrls) {
      dim = new Dimension(totalW, pyrH);
    }
    else {
      dim = new Dimension(totalW, pad);
    }
    transientPanel.setPreferredSize(dim);
    transientPanel.setMinimumSize(dim);
    transientPanel.remove(userDefinedRotPanel);
    transientPanel.repaint();
    revalidate();
    _bigPanel.repaint();
  }

  //...............................

  /**
   *   Makes <em>cursor feedback</em> panels visible.
   */
  protected void addFeedback() {
    feedbackImagePanel.add(feedbackPanel, BorderLayout.NORTH);
    revalidate();
    _bigPanel.repaint();
  }

  /**
   *   Makes <em>cursor feedback</em> panels invisible.
   */
  protected void removeFeedback() {
    feedbackImagePanel.remove(feedbackPanel);
    feedbackImagePanel.repaint();
    revalidate();
    _bigPanel.repaint();
  }

//-------------------------------------------------------------
// methods for opening / closing views
//-------------------------------------------------------------
  /**
   *   Initialises a SectionViewer.
   */
  public void openView() {

    if (_imgV != null) {
      _bigPanel.remove(_imgV);
      _bigPanel.repaint();
      _imageScrollPane.validate();
      _imgV = null;
      _anatBuilder = null;
      _embFileStack = null;
      _xembFileStack = null;
      _embObjStack = null;
      _xembObjStack = null;
      resetFeedbackText();
    }

    _imgV = new WlzImgView();
    _bigPanel.add(_imgV);

    if (_OBJModel != null) {
      _OBJModel = null;
    }
    if (_VSModel != null) {
      _VSModel = null;
    }

    _openingView = true;

    try {
      Method M1 = null;
      M1 = _parent.getClass().getMethod("getOBJModel", null);
      _OBJModel = (WlzObjModel) M1.invoke(_parent, null);
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("getOBJModel: no such method");
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

    _VSModel = new ViewStructModel(_OBJModel);
    _OBJModel.makeSection(_VSModel.getViewStruct());
    setDistLimits(0.0);
    connectAdaptors_1();
    // set WSetters according to which view is set
    resetGUI();
    /*
    try {
      Method M1 = null;
      M1 = _parent.getClass().getMethod("getSVTitleText", null);
      setTitleText( (String) M1.invoke(_parent, null));
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("getSVTitleText: no such method");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }
    //setTitleText(_parent.getTitleText());
    */
    setViewTitle();
    try {
      Method M1 = null;
      M1 = _parent.getClass().getMethod("getAnatomyBuilder", null);
      _anatBuilder = (AnatomyBuilder) M1.invoke(_parent, null);
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("getAnatomyBuilder: no such method 1");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

    Thread t = new Thread(rStack);
    t.start();

    //printModelInfo();

    /* to show intersection of this view on existing views */
    _VSModel.fireChange();

    _openingView = false;

  } // openView()

//---------------------------------------
  /**
   *   Tidies up when a SectionViewer is closed.
   */
  protected void closeView() {
    if (WI2P_1 != null)
      _imgV.removeChangeListener(WI2P_1);
    if (WI2G_1 != null)
      _imgV.removeChangeListener(WI2G_1);
    _imgV = null;
    _embFileStack = null;
    _xembFileStack = null;
    _embObjStack = null;
    _xembObjStack = null;
    _OBJModel = null;
  } // closeView()

//-------------------------------------------------------------
// set max & min for distSetter
// (needs to be done when grey level file is opened &
// after every angle change)
//-------------------------------------------------------------
  /**
   *   Sets <em>max</em> & <em>min</em> limits for <em>distSetter</em>
   *   and sets its current value.
   *   <br>Must be called at initialisation and whenever the 
   *   <em>yaw pitch</em> or <e>roll</em> angles are changed.
   *   @param val the value to set the <em>distSetter</em>.
   */
  protected void setDistLimits(double val) {

    double results[];
    double min = 0.0;
    double max = 0.0;

    boolean state = false;
    boolean sliderState = false;

    results = _OBJModel.getMaxMin(_VSModel.getViewStruct());
    min = results[5];
    max = results[2];

    sliderState = distSetter.isSliderEnabled();
    state = distSetter.isEnabled();

    distSetter.setEnabled(false);
    distSetter.setSliderEnabled(true);
    distSetter.setMin(min);
    distSetter.setMax(max);
    distSetter.setValue(val);
    distSetter.setEnabled(state);
    distSetter.setSliderEnabled(sliderState);
  }

//-------------------------------------------------------------
// set the title text
//-------------------------------------------------------------
/*
  public void setTitleText(String title) {
    titleText.setText(title);
  }
*/

//-------------------------------------------------------------
// re-set Feedback text
//-------------------------------------------------------------
  /**
   *   Clears the text fields associated with <em>cursor feedback</em>.
   *   <ul>
   *   <li>position (x,y,z)</li>
   *   <li>grey value</li>
   *   <li>mouse-click anatomy pathname</li>
   *   </ul>
   */
  protected void resetFeedbackText() {
    resetPosGreyText();
    resetAnatomyText();
  }

//-------------------------------------------------------------
// re-set pos & grey text
//-------------------------------------------------------------
  /**
   *   Clears the feedback text fields associated with
   *   <em>position</em> & <em>grey value</em>.
   */
  protected void resetPosGreyText() {
    xyzTextField.setText("-----");
    valueTextField.setText("-");
  }

//-------------------------------------------------------------
// re-set anatomy text
//-------------------------------------------------------------
  /**
   *   Clears the feedback text fields associated with
   *   <em>mouse-click anatomy</em>.
   */
  protected void resetAnatomyText() {
    anatomyTextField.setText("-----");
  }

//-------------------------------------------------------------
// print model info
//-------------------------------------------------------------
  /**
   *   Debugging utility.
   *   <br>Displays information relating to the underlying models.
   *   <ul>
   *   <li>Woolz Facts</li>
   *   <li>3D Bounding Box</li>
   *   <li>View Structure</li>
   *   </ul>
   */
  protected void printModelInfo() {
    if (_OBJModel != null) {
      _OBJModel.printFacts();
      _OBJModel.printBoundingBox();
      if (_VSModel != null) {
        _VSModel.printViewStructure();
      }
    }
  }

//-------------------------------------------------------------
// print Woolz Facts
//---------------------------------------
  /**
   *   Debugging utility.
   *   <br>Displays information relating to the given Woolz object.
   *   <ul>
   *   <li>Woolz Facts</li>
   *   </ul>
   */
  public void printFacts(WlzObject obj) {
    System.out.println("WLZ FACTS");
    String facts[] = {
        ""};
    try {
      WlzObject.WlzObjectFacts(obj, null, facts, 1);
      System.out.println(facts[0]);
    }
    catch (WlzException e) {
      System.out.println("WlzException #6");
      System.err.println(e);
    }
    System.out.println("*******************************");
  }

//-------------------------------------------------------------
// get anatomy at mouse click
//-------------------------------------------------------------
  /**
   *   Returns the full name to the anatomy component (if any)
   *   at the most recent mouse click / drag.
   *   @param pt3d the position of the most recent mouse click / drag
   *   in the coordinate space of the 3D Woolz object.
   *   @return the full name to the anatomy component, or "".
   */
  protected String getAnatomyAt(double[] pt3d) {

    _maybeFil = new Stack();
    _maybeObj = new Stack();
    File thisFile;
    WlzObject obj = null;
    WlzIBox3 bbox = null;
    int plane = (int) pt3d[2];
    int line = (int) pt3d[1];
    int kol = (int) pt3d[0];
    double smallestVol = 1000000000000000.0;
    String ret = "";
    int i = 0;

    //look through embryo components
    if (_embFileStack != null) {
      int len = _embFileStack.size();
      for (i = 0; i < len; i++) {
        try {
          obj = (WlzObject) _embObjStack.elementAt(i);
          if (WlzObject.WlzInsideDomain(obj, plane, line, kol) != 0) {
            thisFile = (File) _embFileStack.elementAt(i);
            _maybeFil.push(thisFile);
            _maybeObj.push(obj);
          }
        }
        catch (WlzException e) {
          System.out.println("WlzException #1");
          System.out.println(e.getMessage());
        }
      }
    }

    //look through external components
    if (_xembFileStack != null) {
      int xlen = _xembFileStack.size();
      for (i = 0; i < xlen; i++) {
        try {
          obj = (WlzObject) _xembObjStack.elementAt(i);
          if (WlzObject.WlzInsideDomain(obj, plane, line, kol) != 0) {
            thisFile = (File) _xembFileStack.elementAt(i);
            _maybeFil.push(thisFile);
            _maybeObj.push(obj);
          }
        }
        catch (WlzException e) {
          System.out.println("WlzException #2");
          System.out.println(e.getMessage());
        }
      }
    }

    int indx = 0;
    int len = _maybeObj.size();
    //System.out.println(len+" possible anatomy objects");
    double X = 0.0;
    double Y = 0.0;
    double Z = 0.0;
    if (len > 0) {
      // find the smallest wlz object
      for (i = 0; i < len; i++) {
        obj = (WlzObject) _maybeObj.elementAt(i);
        try {
          bbox = WlzObject.WlzBoundingBox3I(obj);
          X = bbox.xMax - bbox.xMin;
          Y = bbox.yMax - bbox.yMin;
          Z = bbox.zMax - bbox.zMin;
        }
        catch (WlzException e) {
          System.out.println("WlzException #?");
          System.out.println(e.getMessage());
        }

        double vol = X * Y * Z;
        if (vol < smallestVol) {
          smallestVol = vol;
          indx = i;
        }
      } // for

      int end = 0;
      StringBuffer temp = new StringBuffer(
          ( (File) _maybeFil.elementAt(indx)).getAbsolutePath());

      // remove /net/... and '.wlz'
      len = AnatomyBuilder._pathLengthToAnatomy;
      temp.replace(0, len, "");
      len = temp.lastIndexOf(SLASH);
      temp.replace(len, temp.length(), "");
      len = temp.lastIndexOf(SLASH);
      if (len != -1) {
        temp.replace(len, temp.length(), capitalise(temp.substring(len)));
        ret = temp.toString();
      }
      else {
        ret = capitalise(temp.toString());
      }
    }

    return ret;
  }

  //-------------------------------------------------------------
  /*
   *   Checks that the given filename ends with ".wlz".
   *   <br>Should also check that the parent directory has the 
   *   same name as the file.
   *   @param dir not used.
   *   @param fil the filename to check.
   *   @return true if the filename ends with ".wlz".
  protected boolean isValidName(String dir, String fil) {
    if (fil.endsWith(".wlz")) {
      return true;
    }
    else {
      return false;
    }
  }
  */

  //-------------------------------------------------------------
  /**
   *   Capitalises the final part of the given filename.
   *   @param name the filename.
   *   @return the filename with the final part capitalised.
   */
  protected String capitalise(String name) {

    String endBit = "";
    StringBuffer buf = new StringBuffer(name);
    String ret = "";

    int lastSlash = -1;

    lastSlash = buf.lastIndexOf(SLASH);
    if (lastSlash != -1) {
      endBit = buf.substring(lastSlash);
      buf.replace(lastSlash, buf.length(), endBit.toUpperCase());
      ret = buf.toString();
    }
    else {
      ret = (buf.toString()).toUpperCase();
    }

    return ret;

  }

//-------------------------------------------------------------
// show anatomy selected from menu
//-------------------------------------------------------------
  /**
   *   Gets an collection of anatomy components from <em>_parent</em> and
   *   calls anatomyFromMenu(Vector vec).
   */
  public void anatomyFromMenu() {

    Vector vec = null;
    if (null != _parent) {
      try {
        Method M1 = null;
        M1 = _parent.getClass().getMethod("getAnatomyElements", null);
        vec = (Vector) M1.invoke(_parent, null);
      }
      catch (InvocationTargetException e) {
        System.out.println(e.getMessage());
      }
      catch (NoSuchMethodException ne) {
        System.out.println("getAnatomyElements: no such method 1");
        System.out.println(ne.getMessage());
      }
      catch (IllegalAccessException ae) {
        System.out.println(ae.getMessage());
      }
      anatomyFromMenu(vec);
    }
  }

  /**
   *   Gets a section through each of the given 3D anatomy components
   *   using the current View Structure
   *   then calls _imgV.setAnatomyObj().
   *   @param vec the collection of 3D anatomy components.
   */
  public void anatomyFromMenu(Vector vec) {

    int num = vec.size();
    if(num != AnatKey.getNRows()) {
       System.out.println("anatomyFromMenu: vec length "+num+
                          " nrows "+AnatKey.getNRows());
    }
    WlzObject obj2D[] = new WlzObject[num];
    boolean viz[] = new boolean[num];
    Color col[] = new Color[num];

    Enumeration els = vec.elements();

    int i = 0;
    AnatomyElement el = null;

    while(els.hasMoreElements()) {
      el = (AnatomyElement)els.nextElement();
      if (el != null) {
        viz[i] = el.isVisible();
	col[i] = AnatKey.getColor(el.getIndx());
        obj2D[i] = _OBJModel.makeSection(
            el.getObj(),
            _VSModel.getViewStruct());
        if (obj2D[i] != null) {
          if (WlzObject.WlzBndObjGetType(obj2D[i]) ==
              WlzObjectType.WLZ_EMPTY_OBJ) {
            obj2D[i] = null;
          }
        }
      }
      else {
        obj2D[i] = null;
        viz[i] = false;
        col[i] = null;
      }
      i++;
    } // while

    try {
      _imgV.clearAnatomy();
      _imgV.setAnatomyObj(obj2D, viz, col);
      _imgV.repaint();
    }
    catch (WlzException err) {
      System.out.println("WlzException #?");
      System.out.println(err.getMessage());
    }

  }

//-------------------------------------------------------------
// get intersection array
//-------------------------------------------------------------
  /**
   *   Returns an array of lines representing the intersection
   *   of this SectionViewer with other open SectionViewers.
   *   @return the intersection lines.
   */
  protected Line2D.Double[] getIntersectionArr() {

    Line2D.Double intersectionArr[] = null;
    SectionViewer SV = null;
    Vector svVec = null;
    try {
      Method M1 = null;
      M1 = _parent.getClass().getMethod("getOpenViews", null);
      svVec = (Vector) M1.invoke(_parent, null);
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("getOpenViews: no such method 1");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

    WlzThreeDViewStruct vs1 = _VSModel.getViewStruct();
    WlzThreeDViewStruct vs2 = null;

    WlzDVertex2 vtx = null;

    int num = svVec.size();

    // angles are in radians
    // theta is the angle of intersection
    double theta = 0.0;
    double tan = 0.0;
    double xmax[] = new double[1];
    double ymax[] = new double[1];
    double xmin[] = new double[1];
    double ymin[] = new double[1];
    double zarr[] = new double[1];

    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
    double xp = 0.0;
    double yp = 0.0;
    double xTotal = 0.0;
    double yTotal = 0.0;

    intersectionArr = new Line2D.Double[num - 1];
    //..............................
    // get max and min values from view struct
    try {
      WlzObject.Wlz3DViewGetMaxvals(vs1, xmax, ymax, zarr);
      WlzObject.Wlz3DViewGetMinvals(vs1, xmin, ymin, zarr);
      xTotal = xmax[0] - xmin[0];
      yTotal = ymax[0] - ymin[0];

    }
    catch (WlzException e0) {
    }
    //..............................
    int j = 0;
    for (int i = 0; i < num; i++) {
      SV = (SectionViewer) svVec.elementAt(i);
      if (!this.equals(SV)) {
         vs2 = SV.getViewStructModel().getViewStruct();
      }
      //..............................
      if (vs2 != null) {
        try {
          // make sure view structs are presented in this order
	  // vs2, vs1
          // or coords will be wrt the wrong view
          vtx = WlzObject.Wlz3DViewGetIntersectionPoint(vs2, vs1);
        }
        catch (WlzException e1) {
	  /*
          Woolz now gives WLZ_ERR_ALG if sections are parallel
	  */
        }

        if (vtx != null) {
	  /*
	   * work-around for bug in Wlz
	   */
	  WlzDVertex2 vtx2 = adjustedIntersectionPoint(vtx);
          /*
          xp = vtx.vtX;
          yp = vtx.vtY;
	  */
          xp = vtx2.vtX;
          yp = vtx2.vtY;
          xp += xTotal/2.0;
          yp += yTotal/2.0;
          // assumes xp lies between 0 and xtotal
          if ( (xp > 0) && (xp < xTotal)) {
            try {
              // make sure view structs are presented in this order
              // or angle will be wrt the wrong view
              theta = WlzObject.Wlz3DViewGetIntersectionAngle(vs2, vs1);
              tan = Math.tan(theta);
            }
            catch (WlzException e2) {
              System.out.println("getIntersectionArr()");
              System.out.println(e2.getMessage());
            }
            if (Math.abs(tan) > 1.0) {
              y0 = 0.0;
              y1 = yTotal;
              x0 = xp + (y0 - yp) / tan;
              x1 = xp + (y1 - yp) / tan;
            }
            else {
              x0 = 0.0;
              x1 = xTotal;
              y0 = yp + (x0 - xp) * tan;
              y1 = yp + (x1 - xp) * tan;
            }
          } else {
              x0 = 0.0;
              y0 = 0.0;
              x1 = 0.0;
              y1 = 0.0;
	  }
          intersectionArr[j] = new Line2D.Double(x0, y0, x1, y1);
	  if(_debug) printIntersection(intersectionArr[j], j);
	  j++;
          vtx2 = null;
	}
      }
      vs2 = null;
      vtx = null;
      theta = 0.0;
      SV = null;
      x0 = 0.0;
      y0 = 0.0;
      x1 = 0.0;
      y1 = 0.0;

    } // for

    return intersectionArr;
  }

//-------------------------------------------------------------
// get interCol array
//-------------------------------------------------------------
  /**
   *   Returns an array representing the colour of each line 
   *   in the <em>intersection line</em> array..
   *   @return array of Color objects.
   */
  protected Color[] getInterColArr() {

    Color colarr[] = null;
    SectionViewer SV = null;
    Vector svVec = null;
    try {
      Method M1 = null;
      M1 = _parent.getClass().getMethod("getOpenViews", null);
      svVec = (Vector) M1.invoke(_parent, null);
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("getOpenViews: no such method 1");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

    int num = svVec.size();
    colarr = new Color[num - 1];
    int j = 0;
    for (int i = 0; i < num; i++) {
      SV = (SectionViewer) svVec.elementAt(i);
      if (!this.equals(SV)) {
         colarr[j] = SV.getSecColorClt().getBackground();
	 j++;
      }
    }

     return colarr;
  }

//-------------------------
  /**
   *   Updates _imgV with
   *   the current intersection line and colour arrays.
   */
  protected void doIntersection() {
    if (_imgV == null) return;
    if (_showIntersection) {
      _imgV.setIntersectionVec(getIntersectionArr());
      _imgV.setInterColVec(getInterColArr());
    }
    _imgV.repaint();
  }

//-------------------------
  /**
   *   Gets the collection of open SectionViewers from _parent
   *   and calls <em>doIntersection()</em> for each one.
   */
  protected void updateIntersections() {

    SectionViewer SV = null;
    Vector svVec = null;
    try {
      Method M1 = null;
      M1 = _parent.getClass().getMethod("getOpenViews", null);
      svVec = (Vector) M1.invoke(_parent, null);
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("getOpenViews: no such method 1");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

    int num = svVec.size();
    for (int i = 0; i < num; i++) {
      SV = (SectionViewer) svVec.elementAt(i);
      SV.doIntersection();
    }
  }
//-------------------------
  /**
   *   Debugging utility.
   *   <br>Prints out the coordinates
   *   of a given <em>intersection line</em> and its
   *   index in the <em>intersection array</em>.
   *   @param line the intersection line.
   *   @param indx array index for this intersection line.
   */
  protected void printIntersection(Line2D.Double line, int indx) {

     String x1Str;
     String y1Str;
     String x2Str;
     String y2Str;

     if(line == null) return;

     x1Str = Double.toString(line.getX1());
     x2Str = Double.toString(line.getX2());
     y1Str = Double.toString(line.getY1());
     y2Str = Double.toString(line.getY2());

     System.out.print("intersection line: "+indx+" "+x1Str+","+y1Str);
     System.out.println(" --- "+x2Str+","+y2Str);

  }

//-------------------------
  /**
   *   Updates _imgV with the coordinates of the fixed point
   *   if it should be visible in the current section.
   */
  protected void doShowFixedPoint() {

    double xa[] = new double[1];
    double ya[] = new double[1];
    double za[] = new double[1];

    double pt3d[] = new double[3];

    int intersect = -1;

    if (_imgV == null)
      return;
    if (!_fixedPoint)
      return;

    _VSModel.getDist(za);
    if (za[0] == 0.0) { // fixed point is in this plane

      _VSModel.getFixedPoint(xa, ya, za);
      pt3d[0] = xa[0];
      pt3d[1] = ya[0];
      pt3d[2] = za[0];
      double vals[] = _OBJModel.get2DPoint(
          pt3d,
          _VSModel.getViewStruct());
      try {
        intersect = WlzObject.WlzInsideDomain(_OBJModel.getSection(),
                                              0.0,
                                              vals[1],
                                              vals[0]);
      }
      catch (WlzException e) {
        System.out.println("doShowFixedPoint");
        System.out.println(e.getMessage());
      }

      if (intersect != 0) { // fixed point is in section
        _imgV.setFixedPointVec(vals);
        vals = null;
        pt3d = null;
        xa = null;
        ya = null;
        za = null;
      }
    }
    else {
      removeFixedPoint();
    }
  }

//-------------------------
  /**
   *   Updates _imgV with the coordinates of the 2nd fixed point
   *   if the fixed line should be visible in the current section.
   */
  protected void doShowFixedLine() {

     double dist[] = new double[1];
     double pt3d[] = new double[3];

     int intersect = -1;

     _drawfixLine = true;

     if (_imgV == null)
        return;
     if (!_axisPoint)
        return;

     _VSModel.getDist(dist);
     /*
      * assumes that fixed line will be set in the same plane
      * as the fixed point
      */
     if(dist[0] == 0.0) {

        pt3d = _VSModel.getAxisPoint();
        double vals[] = _OBJModel.get2DPoint(
              pt3d,
              _VSModel.getViewStruct());
        try {
           intersect = WlzObject.WlzInsideDomain(_OBJModel.getSection(),
                 vals[2],
                 vals[1],
                 vals[0]);
        }
        catch (WlzException e) {
           System.out.println("doShowFixedLine");
           System.out.println(e.getMessage());
        }

        if (intersect != 0) { // axis point is in section
           _imgV.setAxisPointArr(vals);
           vals = null;
           pt3d = null;
        }
     }
     else {
        removeAxisPoint();
     }
  }

//-------------------------------------------------------------
// remove overlay
//-------------------------------------------------------------
  /**
   *   Causes <em>mouse-click anatomy</em> to be cleared from screen.
   */
  protected void removeOverlay() {
    _imgV.clearOverlay();
    _maybeObj = null;
    _maybeFil = null;
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove intersection
//-------------------------------------------------------------
  /**
   *   Causes <em>intersection lines</em> to be cleared from screen.
   */
  protected void removeIntersection() {
    _imgV.clearIntersection();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove threshold
//-------------------------------------------------------------
  /**
   *   Causes the <em>thresholded</em> region to be cleared from screen.
   */
  protected void removeThreshold() {
    _imgV.clearThreshold();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove anatomy
//-------------------------------------------------------------
  /**
   *   Causes <em>anatomy from menu</em> to be cleared from screen.
   */
  public void removeAnatomy() {
    _imgV.clearAnatomy();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove threshConstraint
//-------------------------------------------------------------
  /**
   *   Causes <em>threshold constraint</em> contour to be cleared from screen.
   */
  protected void removeThreshConstraint() {
    _imgV.clearThreshConstraint();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove fixedPoint
//-------------------------------------------------------------
  /**
   *   Causes the <em>fixed point</em> to be cleared from screen.
   */
  protected void removeFixedPoint() {
    _imgV.enableFixedPoint(false);
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove axisPoint
//-------------------------------------------------------------
  /**
   *   Causes the <em>2nd fixed point</em> to be cleared from screen.
   *   <p>Obsolete.
   *   <br>(The 2nd fixed point doesn't get displayed now.)
   */
  protected void removeAxisPoint() {
    _imgV.enableAxisPoint(false);
    _imgV.enableAxis(false);
    _imgV.repaint();
  }

//-------------------------------------------------------------
// set view title
//-------------------------------------------------------------
  /**
   *   Updates the SectionViewer title with
   *   the current values of <em>pitch yaw</em> & <em>roll</em>.
   */
  protected void setViewTitle() {
    // assume all sliders are the same type (double etc)
    int type = pitchSetter.getType();
    Vector pVec = null;
    Vector yVec = null;
    Vector rVec = null;
    String pitch = "";
    String yaw = "";
    String roll = "";
    pVec = pitchSetter.getValue();
    yVec = yawSetter.getValue();
    rVec = rollSetter.getValue();
    switch (type) {
      case INTEGER:
        pitch = ( (Integer) pVec.elementAt(0)).toString();
        yaw = ( (Integer) yVec.elementAt(0)).toString();
        roll = ( (Integer) rVec.elementAt(0)).toString();
        break;
      case FLOAT:
        pitch = ( (Float) pVec.elementAt(0)).toString();
        yaw = ( (Float) yVec.elementAt(0)).toString();
        roll = ( (Float) rVec.elementAt(0)).toString();
        break;
      case DOUBLE:
        pitch = ( (Double) pVec.elementAt(0)).toString();
        yaw = ( (Double) yVec.elementAt(0)).toString();
        roll = ( (Double) rVec.elementAt(0)).toString();
        break;
      default:
        break;
    }
    String viewTitle = trimString(pitch) + " | " +
                       trimString(yaw) + " | " +
                       trimString(roll);
    try {
      Method M1 = null;
      M1 = _frame.getClass().getMethod(
	    "setTitle",
            new Class[]{viewTitle.getClass( )});
      M1.invoke(_frame, new Object[] {viewTitle});
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("setTitle: no such method 2");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

  }

//-------------------------------------------------------------
// remove decimal point and anything past decimal point
//-------------------------------------------------------------
  /**
   *   Utility to remove a decimal point and 
   *   anything following the decimal point
   *   from the given String.
   *   @param str the String to be trimmed.
   *   @return the trimmed String.
   */
  protected String trimString(String str) {

     String ret = "";
     StringBuffer strBuf = new StringBuffer(str);
     int indx = strBuf.indexOf(".");

     if(indx == -1) return str;

     try {
        if(strBuf.length() >= indx+1) {
           ret = strBuf.substring(0, indx+2);
        } else {
           ret = strBuf.substring(0, indx);
        }
     }
     catch(StringIndexOutOfBoundsException e) {
     }

     return ret;
  }

//-------------------------------------------------------------
// save image as jpeg
//-------------------------------------------------------------
  /**
   *   Saves the composite screen image to a jpg file.
   *   @param str filename for the jpg image.
   */
  protected void saveImage(String str) {
    File fil = new File(str);
    BufferedImage img = null;
    /* boolean needed because 3D apps have code in WlzImgView */
    boolean showIntersectionLines = true;

    img = _imgV.getComponentBufferedImage(showIntersectionLines);

    try {
      FileOutputStream fout = new FileOutputStream(fil);
      JPEGImageEncoder JIE = JPEGCodec.createJPEGEncoder(fout);
      JIE.encode(img);
    }
    catch (FileNotFoundException e1) {
      System.out.println(e1.getMessage());
    }
    catch (ImageFormatException e2) {
      System.out.println(e2.getMessage());
    }
    catch (IOException e3) {
      System.out.println(e3.getMessage());
    }
  }

//-------------------------------------------------------------
// save view settings
//-------------------------------------------------------------
  /**
   *   Saves the current <em>ViewStruct</em> and zoom setting to an xml file.
   *   <br>The xml file may be read back in using <em>loadViewFromXML</em>
   *   @param str filename for the xml data.
   */
  protected void saveViewAsXML(String fil) {

    String tab = "   ";
    String str = "";
    double val1[] = new double[1];
    double val2[] = new double[1];
    double val3[] = new double[1];
    double rtod = 180.0 / Math.PI;

    try {
      PrintWriter pout = new PrintWriter(
          new FileWriter(fil));

      pout.println("<?xml version=\"1.0\"?>");
      pout.println("<view>");
      //.......................
      pout.print("<zoom>");
      pout.print(Integer.toString(zoomSetter.getValue()));
      pout.println("</zoom>");
      //.......................
      pout.println("<viewStruct>");
      //.......................
      val1[0] = 0.0;
      _VSModel.getDist(val1);
      pout.print(tab);
      pout.print("<dist>");
      str = Double.toString(val1[0]);
      pout.print(str);
      pout.println("</dist>");
      //.......................
      val1[0] = 0.0;
      _VSModel.getTheta(val1);
      pout.print(tab);
      pout.print("<theta>");
      str = Double.toString(val1[0] * rtod);
      pout.print(str);
      pout.println("</theta>");
      //.......................
      val1[0] = 0.0;
      _VSModel.getPhi(val1);
      pout.print(tab);
      pout.print("<phi>");
      str = Double.toString(val1[0] * rtod);
      pout.print(str);
      pout.println("</phi>");
      //.......................
      val1[0] = 0.0;
      _VSModel.getZeta(val1);
      pout.print(tab);
      pout.print("<zeta>");
      str = Double.toString(val1[0] * rtod);
      pout.print(str);
      pout.println("</zeta>");
      //.......................
      val1[0] = 0.0;
      val2[0] = 0.0;
      val3[0] = 0.0;
      _VSModel.getFixedPoint(val1, val2, val3);
      pout.print(tab);
      pout.println("<fixed>");
      pout.print(tab + tab);
      pout.print("<x>");
      str = Double.toString(val1[0]);
      pout.print(str);
      pout.println("</x>");
      pout.print(tab + tab);
      pout.print("<y>");
      str = Double.toString(val2[0]);
      pout.print(str);
      pout.println("</y>");
      pout.print(tab + tab);
      pout.print("<z>");
      str = Double.toString(val3[0]);
      pout.print(str);
      pout.println("</z>");
      pout.print(tab);
      pout.println("</fixed>");
      //.......................
      String str1 = _VSModel.getViewMode();
      str = str1.substring(0, str1.length() - 1);
      pout.print(tab);
      pout.print("<mode>");
      if (str.equals("WLZ_ZETA_MODE") == true) {
        pout.print("absolute");
      }
      else {
        pout.print("up_is_up");
      }
      pout.println("</mode>");
      //.......................
      pout.println("</viewStruct>");
      pout.println("</view>");

      if (pout.checkError() == true) {
        System.out.println("error writing to xml file ");
      }
    }
    catch (FileNotFoundException e1) {
      System.out.println(e1.getMessage());
    }
    catch (IOException e3) {
      System.out.println(e3.getMessage());
    }
  }

//-------------------------------------------------------------
// load view settings
//-------------------------------------------------------------
  /**
   *   Loads <em>ViewStruct</em> and zoom setting from an xml file.
   *   <br>This causes a previously saved section to be displayed.
   *   @param str filename for the xml data.
   */
  protected void loadViewFromXML(File fil) {

    ParseXML parser = new ParseXML(this);
    try {
      parser.doParse(fil);
    }
    catch (Exception e) {
      System.out.println(e.getMessage());
    }

  }

//-------------------------------------------------------------
// convert from (slider) zoom factor to mag
//-------------------------------------------------------------
  /**
   *   Converts from <em>zoom factor</em> to
   *   <em>% magnification</em>.
   *   <br>Zoom factor is displayed on the zoom control whereas
   *   % magnification is used by the program.
   *   @param zf the current zoom factor.
   *   @return % magnification.
   */
  protected double zf2mag(double zf) {
    return zf/100.0;
  }

//-------------------------------------------------------------
// checkFilename for extension
//-------------------------------------------------------------
  /**
   *   Adds the given extension to a filename if it 
   *   doesn't already have it.
   *   <br>Used for jpg and xml files.
   *   @param str the filename.
   *   @param ext the required filename extension,
   *   (not including .).
   *   @return the filename with extension.
   */
  protected String checkFilename(String str, String ext) {

    StringBuffer strbuf;
    String lcaseStr = ext.toLowerCase();
    String ucaseStr = ext.toUpperCase();
    int len = 0;

    if ( (str.endsWith(lcaseStr)) || (str.endsWith(ucaseStr))) {
      return str;
    }
    else {
      strbuf = new StringBuffer(str);
      strbuf.append("." +lcaseStr);
      return strbuf.toString();
    }
  }

//-------------------------------------------------------------
// set view mode
//-------------------------------------------------------------
  /**
   *   Sets the <em>View Structure</em> with the given
   *   view mode and initialises the rotation sliders accordingly.
   *   @param mode the view mode to set.
   */
  public void setViewMode(String mode) {
     if (_VSModel != null) {

        double[] zeta = null;
        double roll = 0.0;

        if (mode.equals("up is up")) {
           _fixedLineRotation = false;
           resetFixedLine(false);
	   _VSModel.setViewMode(mode);
           controlMenu_2_1.setSelected(true);
	   enableFPMenus(true);
           yawSetter.setSliderEnabled(true);
           yawSetter.setEnabled(true);
           pitchSetter.setSliderEnabled(true);
           pitchSetter.setEnabled(true);
           rollSetter.setSliderEnabled(true);
	   zeta = new double[1];
	   _VSModel.getZeta(zeta);
	   roll = zeta[0] * 180.0 / Math.PI;
           rollSetter.setValue(roll);
           rollSetter.setSliderEnabled(false);
           rotSetter.setEnabled(false);
           rotSetter.setSliderEnabled(true);
           rotSetter.setValue(0.0);
           rotSetter.setSliderEnabled(false);
        }
        else if (mode.equals("absolute")) {
           _fixedLineRotation = false;
           resetFixedLine(false);
	   _VSModel.setViewMode(mode);
           controlMenu_2_2.setSelected(true);
	   enableFPMenus(true);
           yawSetter.setSliderEnabled(true);
           yawSetter.setEnabled(true);
           pitchSetter.setSliderEnabled(true);
           pitchSetter.setEnabled(true);
           rollSetter.setSliderEnabled(true);
           rollSetter.setEnabled(true);
           rotSetter.setEnabled(false);
           rotSetter.setSliderEnabled(true);
           rotSetter.setValue(0.0);
           rotSetter.setSliderEnabled(false);
        }
        else if(mode.equals("fixed line")) {
           _fixedLineRotation = true;
           controlMenu_2_3.setEnabled(true);
           controlMenu_2_3.setSelected(true);
           controlMenu_2_3.setEnabled(false);
	   enableFPMenus(false);
           yawSetter.setSliderEnabled(false);
           yawSetter.setEnabled(false);
           pitchSetter.setSliderEnabled(false);
           pitchSetter.setEnabled(false);
           rollSetter.setSliderEnabled(false);
           rollSetter.setEnabled(false);
           rotSetter.setSliderEnabled(true);
           rotSetter.setEnabled(true);
	   _VSModel.setViewMode(mode);
        }
     }
  }

//-------------------------------------------------------------
// enable or disable fixed point menus depending upon view mode
//-------------------------------------------------------------
  /**
   *   Enables or disables <em>Fixed Point</em> menus 
   *   depending upon the given state.
   *   @param state true if menus are to be enabled.
   */
  private void enableFPMenus(boolean state) {
     controlMenu_3_1.setEnabled(state);
     controlMenu_3_2.setEnabled(state);
     controlMenu_3_3.setEnabled(state);
     controlMenu_4_1.setEnabled(state);
     controlMenu_4_2.setEnabled(state);
     controlMenu_4_3.setEnabled(state);
  }

//-------------------------------------------------------------
// get view structure model etc
//-------------------------------------------------------------
  /**
   *   Returns the current <em>ViewStructModel</em>.
   *   @return _VSModel.
   */
  public ViewStructModel getViewStructModel() {
    return _VSModel;
  }

  /**
   *   Returns the current <em>WlzObjModel</em>.
   *   @return _OBJModel.
   */
  protected WlzObjModel getWlzObjModel() {
    return _OBJModel;
  }

  /**
   *   Returns the current <em>WlzImgView</em>.
   *   @return _imgV.
   */
  protected WlzImgView getImageView() {
    return _imgV;
  }

  /**
   *   Returns the <em>zoom</em> control.
   *   @return zoomSetter.
   */
  protected Zoom getZoomSetter() {
    return zoomSetter;
  }

  /**
   *   Returns the <em>dist</em> control.
   *   @return distSetter.
   */
  protected WSetter getDistSetter() {
    return distSetter;
  }

  /**
   *   Returns the <em>yaw</em> control.
   *   @return yawSetter.
   */
  protected WSetter getYawSetter() {
    return yawSetter;
  }

  /**
   *   Returns the <em>pitch</em> control.
   *   @return pitchSetter.
   */
  protected WSetter getPitchSetter() {
    return pitchSetter;
  }

  /**
   *   Returns the <em>roll</em> control.
   *   @return rollSetter.
   */
  protected WSetter getRollSetter() {
    return rollSetter;
  }

  /**
   *   Returns the <em>colour chooser</em> control.
   *   @return secColorClt.
   */
  protected JButton getSecColorClt() {
    return secColorClt;
  }

//-------------------------------------------------------------
// enable / disable thresholding & threshConstraint
//-------------------------------------------------------------
  /**
   *   Toggles <em>threshold</em> mode on or off.
   *   <br>When in thresholdmode, pressing the mouse button
   *   over the grey level image will use the current grey value
   *   as the mid-point of the threshold range.
   *   <br>Subsequent dragging of the mouse to the right or left will
   *   increase or decrease the threshold window respectively.
   *   <p>Not in use.
   *   @param state true if thresholding is enabled.
   */
  protected void enableThresholding(boolean state) {
    _thresholding = state;
  }

//-------------------------------------------------------------
  /**
   *   Toggles <em>threshold constraint</em> mode on or off.
   *   <br>When in threshold constraint mode, pressing and dragging the mouse
   *   over the grey-level image will draw a contour.
   *   When the mouse button is released the contour automatically closes.
   *   <br>The constraint region is the inside of the closed contour.
   *   <p>Not in use.
   *   @param state true if threshold constraint is enabled.
   */
  protected void enableThreshConstraint(boolean state) {
    _threshConstraint = state;
    if (state) {
      if (!thresholdMenu_3state) {
        thresholdMenu_3.doClick();
      }
    }
  }

//-------------------------------------------------------------
  /**
   *   Disables <em>mouse-click anatomy</em> mode.
   */
  protected void disableAnatomy() {
    _anatomy = false;
  }

  /**
   *   Enables <em>mouse-click anatomy</em> mode.
   */
  protected void enableAnatomy() {
    _anatomy = true;
  }

//-------------------------------------------------------------
// set anatomy text field for menu selected anatomy
//-------------------------------------------------------------
  /**
   *   Sets the <em>feedback anatomy text</em> to the current 
   *   <em>mouse-click anatomy</em> or "-----" if none.
   *   <br>Called from SVParent2D when <em>Anatomy Key</em> is cleared..
   *   @param str the String to set.
   *   <p>Note: It might be better to get rid of str as SVParent2D just
   *   sets it to SectionViewer._anatomyStr.
   */
  public void setAnatomyText(String str) {
    resetAnatomyText();
    anatomyTextField.setText(str);
  }

//-------------------------------------------------------------
// reset fixed point
//-------------------------------------------------------------
  /**
   *   Sets <em>Fixed Point</em> to its default value and
   *   displays it if showFP is true.
   *   @param showFP true if fixed point is to be displayed.
   */
  public void resetFixedPoint(boolean showFP) {

    double vals[] = null;

    // reset fixed point
    vals = _VSModel.getInitialFixedPoint();
    _VSModel.setFixedPoint(vals[0], vals[1], vals[2]);
    vals = null;
    setDistLimits(0.0);
    if (showFP) {
      Runnable fpClick = new Runnable() {
        public void run() {
          if (!showMenu_4state)
            showMenu_4.doClick();
        }
      };
      SwingUtilities.invokeLater(fpClick);
    }
  }

//-------------------------------------------------------------
// reset fixed line
//-------------------------------------------------------------
  /**
   *   Sets <em>2nd Fixed Point</em> to the original 
   *   <em>Fixed Point</em> and
   *   displays the fixed line if showAP is true.
   *   @param showAP true if fixed line is to be displayed.
   */
  public void resetFixedLine(boolean showAP) {

     double vals[] = null;
     final boolean show = showAP; // required for use in inner class

     // reset fixed line
     vals = _VSModel.getFixedPoint();
     _VSModel.setAxisPoint(vals[0],
	   vals[1],
	   vals[2]);
     vals = null;
     Runnable apClick = new Runnable() {
	public void run() {
	   if ((!showMenu_5state && show) ||
		 (showMenu_5state && !show)) {
	      showMenu_5.doClick();
	   }
	}
     };
     SwingUtilities.invokeLater(apClick);
  }

//-------------------------------------------------------------
// reset user interface controls
//-------------------------------------------------------------
  /**
   *   Sets the SectionViewer controls to default values.
   */
  protected void resetGUI() {
    if (_debug)
      System.out.println("enter resetGUI");

    setSVEnabled(false);
    zoomSetter.setEnabled(true);
    distSetter.setEnabled(true);
    pitchSetter.setEnabled(true);
    yawSetter.setEnabled(true);
    rollSetter.setEnabled(true);
    rotSetter.setEnabled(true);

    if (viewType.equals("XY")) {
      pitchSetter.setValue(0.0);
      yawSetter.setValue(0.0);
      rollSetter.setValue(0.0);
    }
    else if (viewType.equals("YZ")) {
      pitchSetter.setValue(90.0);
      yawSetter.setValue(0.0);
      rollSetter.setValue(90.0);
    }
    else if (viewType.equals("ZX")) {
      pitchSetter.setValue(90.0);
      yawSetter.setValue(90.0);
      rollSetter.setValue(90.0);
    }
    else {
      System.out.println("unknown view type in resetGUI(): " + viewType);
    }
    distSetter.setValue(0.0);
    zoomSetter.setValue(100);
    rotSetter.setValue(0.0);
    distSetter.setSliderEnabled(true);
    pitchSetter.setSliderEnabled(true);
    yawSetter.setSliderEnabled(true);
    rollSetter.setSliderEnabled(false);
    rotSetter.setSliderEnabled(false);
    setViewMode("up is up");
    removeThreshold();
    removeThreshConstraint();
    enableThresholding(false);
    enableThreshConstraint(false);
    disableAnatomy();
    resetFeedbackText();
    resetFixedPoint(false);
    resetFixedLine(false);
    setSVEnabled(true);

    if (_debug)
      System.out.println("exit resetGUI");
  }

//-------------------------------------------------------------
// constraint to Wlz
//-------------------------------------------------------------
  /**
   *   Converts the region defined by a <em>threshold constraint</em>
   *   to a Woolz object.
   *   <br>This is required for constrained thresholding.
   *   @param bbox2 the 2D bounding box of the current section.
   */
  protected void constraintToWlz(WlzIBox2 bbox2) {

    GeneralPath GP = null;
    PathIterator PI = null;
    WlzIVertex2 vertArray[] = null;

    int WLZ_POLYGON_INT = 1;
    double ofsX = (bbox2.xMax - bbox2.xMin) / 2.0;
    double ofsY = (bbox2.yMax - bbox2.yMin) / 2.0;

    GP = _imgV.getThreshConstraint();
    PI = GP.getPathIterator(null);
    int tot = 0;
    while (!PI.isDone()) {
      PI.next();
      tot++;
    }
    int maxv = tot - 1;
    PI = null;
    PI = GP.getPathIterator(null);
    float segment[] = new float[6];

    /* the first segment has coordinates 0,0 and is ignored */
    vertArray = new WlzIVertex2[maxv];
    int j = 0;
    for (int i = 0; i < tot; i++) {
      if ( (PI.currentSegment(segment) == PI.SEG_LINETO) ||
          (PI.currentSegment(segment) == PI.SEG_CLOSE)) {

        vertArray[j] = new WlzIVertex2(
            (int) (segment[0] - ofsX),
            (int) (segment[1] - ofsY));
        j++;
      }
      PI.next();
    }
    if (!PI.isDone()) {
      System.out.println("iterator not done");
      System.out.println("array length " + vertArray.length);
    }

    int copy = 1; // allocate space and copy
    int fillType = 0; // WLZ_SIMPLE_FILL
    WlzPolygonDomain WPD = null;
    int transformType = 2; // WLZ_TRANSFORM_2D_TRANS

    try {
      WPD = WlzObject.WlzMakePolygonDomain(
          WLZ_POLYGON_INT,
          maxv,
          vertArray,
          maxv,
          copy);

      _constraintObj = WlzObject.WlzPolyToObj(WPD, fillType);

    }
    catch (WlzException we) {
      System.out.println("WlzMakePolygonDomain");
      System.out.println(we.getMessage());
    }
  }

//-------------------------------------------------------------
// close
//-------------------------------------------------------------
  /**
   *   closes a SectionViewer.
   */
  protected void close() {
    try {
      Method M1 = null;
      M1 = _parent.getClass().getMethod(
	    "removeView",
            new Class[]{this.getClass( )});
      M1.invoke(_parent, new Object[] {this});
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("removeView: no such method 2");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

    closeView();
    updateIntersections();

    try {
      Method M1 = null;
      M1 = _frame.getClass().getMethod(
	    "dispose", null);
      M1.invoke(_frame, null);
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("dispose: no such method ");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }
  }

  /**
   *   Work-around for bug in Wlz.
   *   <br>Wlz3DViewGetIntersectionPoint gives a shifted answer
   *   if the fixed point is changed. This causes intersection lines
   *   to jump in the view with the changed fixed point
   *   @param vtx the intersection point returned by Woolz.
   *   @return the modified intersection point.
   */
  private WlzDVertex2 adjustedIntersectionPoint(WlzDVertex2 vtx) {

     if(_debug) System.out.println("entering adjustIntersectionPoint");

     double[] fpInitial = null;
     double[] fpNew = new double[3];
     double[] x = new double[1];
     double[] y = new double[1];
     double[] z = new double[1];
     double[] fp2DOrig = null;
     double[] fp2DNew = null;

     WlzDVertex2 ret = vtx;
     WlzThreeDViewStruct vs = _VSModel.getViewStruct();

     _VSModel.getFixedPoint(x, y, z);
     fpNew[0] = x[0];
     fpNew[1] = y[0];
     fpNew[2] = z[0];
     fpInitial = _VSModel.getInitialFixedPoint();

     fp2DOrig = _OBJModel.get2DPoint(fpInitial, vs);
     fp2DNew = _OBJModel.get2DPoint(fpNew, vs);

     ret.vtX -= fp2DOrig[0] - fp2DNew[0];
     ret.vtY -= fp2DOrig[1] - fp2DNew[2];

     if(_debug) System.out.println("exiting adjustIntersectionPoint");

     return ret;
  }
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//-------------------------------------------------------------
// connect up the adaptors
//-------------------------------------------------------------
  /** Event handler for the <em>zoom</em> control. */
  ZoomToZoomAdaptor Z2Z_1;

  /** Event handler for the <em>dist</em> control. */
  WSetterToDistAdaptor W2D_1;

  /** Event handler for the <em>pitch</em> control. */
  WSetterToPitchAdaptor W2P_1;

  /** Event handler for the <em>yaw</em> control. */
  WSetterToYawAdaptor W2Y_1;

  /** Event handler for the <em>roll</em> control. */
  WSetterToRollAdaptor W2R_1;

  /** Event handler for the <em>fixed line</em> control. */
  WSetterToFixedRotAdaptor W2FR_1;

  /**
   *   Event handler which changes the section
   *   when the <em>View Structure</em> changes.
   */
  ViewStructModelToViewAdaptor WOM2V_1;

  /**
   *   Event handler which changes
   *   <em>position</em> feedback
   *   when <em>_imgV</em> fires an event.
   */
  WlzImgViewToPosAdaptor WI2P_1;

  /**
   *   Event handler which changes
   *   <em>grey value</em> feedback
   *   when <em>_imgV</em> fires an event.
   */
  WlzImgViewToGreyAdaptor WI2G_1;

  /**
   *   Event handler which changes
   *   <em>mouse-click anatomy</em> feedback
   *   when <em>_imgV</em> fires an event.
   */
  WlzImgViewToAnatAdaptor WI2A_1;

  /**
   *   Event handler that updates <em>_imgV</em> with
   *   the current cursor position in the 2D grey-level image,
   *   then causes <em>_imgV</em> to fire a changeEvent.
   *   <br>Needs to be <em>public</em> as external applications use it.
   */
  public MouseToFBAdaptor M2FB_1; 

  /**
   *   Event handler which changes
   *   <em>fixed point</em> according to cursor position
   *   in grey-level image.
   */
  MouseToFPAdaptor M2FP_1;

  /**
   *   Event handler which generates a
   *   <em>threshold</em> region according to
   *   cursor position & movement in grey-level image.
   */
  MouseToThresholdAdaptor M2T_1;

  /**
   *   Event handler which generates a
   *   <em>threshold constraint</em> contour according to
   *   cursor position & movement in grey-level image.
   */
  MouseToThreshConstraintAdaptor M2TC_1;

  /**
   *   Event handler which fires an event if 
   *   the mouse is clicked in a SectionViewer.
   */
  FocusAdaptor FA_1;

  // adaptors for grey level
  /**
   *   Creates new event handlers and registers them with the 
   *   appropriate event-firing objects.
   *   <br>Any existing event handlers are removed prior to registering
   *   the new ones.
   */
  protected void connectAdaptors_1() {
    if (_debug)
      System.out.println("enter connectAdaptors_1");

    setSVEnabled(false);
    //...............................
    // WSetter to model adaptors
    //...............................
    if (Z2Z_1 != null)
      zoomSetter.removeChangeListener(Z2Z_1);
    Z2Z_1 = new ZoomToZoomAdaptor(zoomSetter, _imgV);
    zoomSetter.addChangeListener(Z2Z_1);
    //...............................
    if (distSetter == null)
      System.out.println("distSetter = null");
    if (_VSModel == null)
      System.out.println("_VSModel = null");
    if (_OBJModel == null)
      System.out.println("_OBJModel = null");
    if (W2D_1 != null)
      distSetter.removeChangeListener(W2D_1);
    W2D_1 = new WSetterToDistAdaptor(distSetter, _VSModel, _OBJModel);
    if (W2D_1 == null)
      System.out.println("W2D_1 = null");
    distSetter.addChangeListener(W2D_1);
    //...............................
    if (W2P_1 != null)
      pitchSetter.removeChangeListener(W2P_1);
    W2P_1 = new WSetterToPitchAdaptor(pitchSetter, _VSModel, _OBJModel);
    pitchSetter.addChangeListener(W2P_1);
    //...............................
    if (W2Y_1 != null)
      yawSetter.removeChangeListener(W2Y_1);
    W2Y_1 = new WSetterToYawAdaptor(yawSetter, _VSModel, _OBJModel);
    yawSetter.addChangeListener(W2Y_1);
    //...............................
    if (W2R_1 != null)
      rollSetter.removeChangeListener(W2R_1);
    W2R_1 = new WSetterToRollAdaptor(rollSetter, _VSModel, _OBJModel);
    rollSetter.addChangeListener(W2R_1);
    //...............................
    if (W2FR_1 != null)
      rotSetter.removeChangeListener(W2FR_1);
    W2FR_1 = new WSetterToFixedRotAdaptor(rotSetter, _VSModel, _OBJModel);
    rotSetter.addChangeListener(W2FR_1);
    //...............................
    // View adaptors
    //...............................
    if (WOM2V_1 != null)
      _VSModel.removeChangeListener(WOM2V_1);
    WOM2V_1 = new ViewStructModelToViewAdaptor(_VSModel, _OBJModel, _imgV);

    _VSModel.addChangeListener(WOM2V_1);
    //...............................
    if (WI2P_1 != null)
      _imgV.removeChangeListener(WI2P_1);
    WI2P_1 = new WlzImgViewToPosAdaptor(_imgV, _OBJModel, _VSModel,
                                        xyzTextField);
    _imgV.addChangeListener(WI2P_1);
    //...............................
    if (WI2G_1 != null)
      _imgV.removeChangeListener(WI2G_1);
    WI2G_1 = new WlzImgViewToGreyAdaptor(_imgV, valueTextField);
    _imgV.addChangeListener(WI2G_1);
    //...............................
    if (WI2A_1 != null)
      _imgV.removeChangeListener(WI2A_1);
    WI2A_1 = new WlzImgViewToAnatAdaptor(_VSModel, _OBJModel, _imgV,
                                         anatomyTextField);
    _imgV.addChangeListener(WI2A_1);
    //...............................
    if (M2FB_1 != null) {
      _imgV.removeMouseListener(M2FB_1);
      _imgV.removeMouseMotionListener(M2FB_1);
    }
    M2FB_1 = new MouseToFBAdaptor(_imgV);
    _imgV.addMouseListener(M2FB_1);
    _imgV.addMouseMotionListener(M2FB_1);
    //...............................
    if (M2FP_1 != null) {
      _imgV.removeMouseListener(M2FP_1);
      _imgV.removeMouseMotionListener(M2FP_1);
    }
    M2FP_1 = new MouseToFPAdaptor(_imgV, _OBJModel, _VSModel);
    _imgV.addMouseListener(M2FP_1);
    _imgV.addMouseMotionListener(M2FP_1);
    //...............................
    if (M2T_1 != null) {
      _imgV.removeMouseListener(M2T_1);
      _imgV.removeMouseMotionListener(M2T_1);
    }
    M2T_1 = new MouseToThresholdAdaptor(_imgV, _OBJModel);
    _imgV.addMouseListener(M2T_1);
    _imgV.addMouseMotionListener(M2T_1);
    //...............................
    if (M2TC_1 != null) {
      _imgV.removeMouseListener(M2TC_1);
      _imgV.removeMouseMotionListener(M2TC_1);
    }
    M2TC_1 = new MouseToThreshConstraintAdaptor(_imgV,
                                                _OBJModel);
    _imgV.addMouseListener(M2TC_1);
    _imgV.addMouseMotionListener(M2TC_1);
    //...............................
    // Focus adaptor
    //...............................
    if (FA_1 != null)
      //this.getContentPane().removeMouseListener(FA_1);
      _imgV.removeMouseListener(FA_1);
      FA_1 = new FocusAdaptor();

      _imgV.addMouseListener(FA_1);
    //...............................
    //--------------------------------------

    setSVEnabled(true);
    if (_debug)
      System.out.println("exit connectAdaptors_1");
  } //connectAdaptors_1()

//-------------------------------------------------------------
  /** 
   *   Toggles the SectionViewer's ability to fire events.
   *   @param state true if the SectionViewer can fire events.
   */
  protected void setSVEnabled(boolean state) {
    _enabled = state;
  }

//-------------------------------------------------------------
  /** 
   *   Tests the SectionViewer's ability to fire events.
   *   @return true if the SectionViewer can fire events.
   */
  protected boolean isSVEnabled() {
    return _enabled;
  }

//-------------------------------------------------------------
  /**
   *   Calculates various parameters required for
   *   <em>fixed line view mode</em> and sets them in _VSModel.
   *   <ul>
   *      <li>fixed line angle</li>
   *      <li>_norm2</li>
   *      <li>_norm3</li>
   *   </ul>
   */
  protected void setFixedLineParams() {
    /*
     * get the coords of fixed point and fixed line in 2D
     * calculate fixed line angle in 2D and set View Struct.
     * calculate vectors norm2 and norm3
     * (see viewFixedPointUtils.c fixed_2_cb())
     */

     double fp3D[];
     double ap3D[];
     double fp2D[];
     double ap2D[];

     double arg = 0.0;
     double angle = 0.0;

     double x = 0.0;
     double y = 0.0;
     double z = 0.0;
     double s = 0.0;

     WlzThreeDViewStruct VS = _VSModel.getViewStruct();
     WlzAffineTransform trans = null;

     double transMatrix[][][] = new double[1][][];
     transMatrix[0] = null;

     fp3D = _VSModel.getFixedPoint();
     ap3D = _VSModel.getAxisPoint();

     /*
     System.out.println("fixed point 3D = ");
     printDoubleArray(fp3D);
     System.out.println("axis point 3D = ");
     printDoubleArray(ap3D);
     */

     fp2D = _OBJModel.get2DPoint(fp3D, VS);
     ap2D = _OBJModel.get2DPoint(ap3D, VS);

     /*
     System.out.println("fixed point 2D = ");
     printDoubleArray(fp2D);
     System.out.println("axis point 2D = ");
     printDoubleArray(ap2D);
     */

     angle = Math.atan2(ap2D[1]-fp2D[1], ap2D[0]-fp2D[0]);
     //System.out.println("angle (degs) = "+angle*180.0/Math.PI);

     _VSModel.setFixedLineAngle(angle);

     /* define the 2 orthogonal axes
        first the normalised direction vector between the fixed points */
     x = ap3D[0] - fp3D[0];
     y = ap3D[1] - fp3D[1];
     z = ap3D[2] - fp3D[2];
     s = Math.sqrt((double) x*x + y*y + z*z + 0.00001);
     if( s <= .001 ){
       x = 1.0;
       y = 0.0;
       z = 0.0;
     } else {
       x /= s; y /= s; z /= s;
     }

     trans = _VSModel.getInverseTransform();

     transMatrix = _VSModel.getTransMatrix(trans);
     /*
     if(transMatrix[0] != null){
        System.out.println(" got transform matrix");
	//printMatrix(transMatrix[0]);
     } else {
        System.out.println("transMatrix[0] == null");
     }
     */

     double mat02 = transMatrix[0][0][2];
     double mat12 = transMatrix[0][1][2];
     double mat22 = transMatrix[0][2][2];

     _VSModel.setNorm2(mat02, mat12, mat22);
     _VSModel.setNorm3(y*mat22 - z*mat12,
                       z*mat02 - x*mat22,
		       x*mat12 - y*mat02);

     /*
     System.out.println("SectionViewer.setFixedLineParams:");
     System.out.println("norm2");
     printDoubleArray(_VSModel.getNorm2());
     System.out.println("norm3");
     printDoubleArray(_VSModel.getNorm3());
     */

  }

//-------------------------------------------------------------
  /**
   *   Debugging utility.
   *   <br>Prints a 3 element array representing a point in 3D space.
   *   @param pnt the 3 element array.
   */
  private void printDoubleArray(double[] pnt) {
     System.out.println("x = "+Double.toString(pnt[0]));
     System.out.println("y = "+Double.toString(pnt[1]));
     System.out.println("z = "+Double.toString(pnt[2]));
     System.out.println("---------------");
  }
//-------------------------------------------------------------
  /**
   *   Debugging utility.
   *   <br>Prints a 3x3 matrix representing rotation.
   *   @param mtrx a 3x3 matrix.
   */
  private void printMatrix(double[][] mtrx) {
     System.out.println("00 = "+Double.toString(mtrx[0][0]));
     System.out.println("01 = "+Double.toString(mtrx[0][1]));
     System.out.println("02 = "+Double.toString(mtrx[0][2]));
     System.out.println("10 = "+Double.toString(mtrx[1][0]));
     System.out.println("11 = "+Double.toString(mtrx[1][1]));
     System.out.println("12 = "+Double.toString(mtrx[1][2]));
     System.out.println("20 = "+Double.toString(mtrx[2][0]));
     System.out.println("21 = "+Double.toString(mtrx[2][1]));
     System.out.println("22 = "+Double.toString(mtrx[2][2]));
     System.out.println("---------------");
  }

//-------------------------------------------------------------
  /**
   *   Constrains <em>Fixed Point</em> and <em>2nd Fixed Point</em>
   *   to be within the 3D bounding box.
   *   @param pnt a 3 element array representing the position of the cursor
   *   when setting <em>Fixed Point</em> or <em>2nd Fixed Point</em>.
   */
  private void clipPoint3D(double[] pnt) {

     WlzIBox3 bBox = null;
     bBox = _OBJModel.getBBox();

     if(bBox == null) return;

     pnt[0] = pnt[0] < bBox.xMin ? bBox.xMin : pnt[0];
     pnt[0] = pnt[0] > bBox.xMax ? bBox.xMax : pnt[0];
     pnt[1] = pnt[1] < bBox.yMin ? bBox.yMin : pnt[1];
     pnt[1] = pnt[1] > bBox.yMax ? bBox.yMax : pnt[1];
     pnt[2] = pnt[2] < bBox.zMin ? bBox.zMin : pnt[2];
     pnt[2] = pnt[2] > bBox.zMax ? bBox.zMax : pnt[2];

  }

//-------------------------------------------------------------
  /* checks that dist = 0.0 when setting fixed line */
  /**
   *   True if the <em>Fixed Point</em> and <em>2nd Fixed Point</em>
   *   are in the same plane when <em>2nd Fixed Point</em>
   *   is being changed.
   *   @return true or false.
   */
  private boolean checkDist() {
     boolean ret = false;
     double[] dArr = null;
     String message = "";
     int reply = 0;

     dArr = new double[1];

     _VSModel.getDist(dArr);
     if(dArr[0] == 0.0) {
        ret = true;
     } else {
        message = "the fixed line end-point must be in the same plane as the fixed point (with dist = 0.0)\nIf you want the fixed line to be in this plane please change the fixed point first.";
	JOptionPane.showMessageDialog(
	     null,
	     message,
	     "setting fixed line",
	     JOptionPane.INFORMATION_MESSAGE);
     }

     return ret;
  }

//-------------------------------------------------------------
  /**
   *   Inverts the grey values of the grey-level image.
   *   Code was moved from <em>invertButton</em> handler
   *   so that applications can override this method.
   */
  public void invertSection() {
    // reverses the grey-level image
    if(!_inverted) {
       invertButton.setBackground(Color.black);
       invertButton.setForeground(Color.white);
       _imgV.setInverted(true);
    } else {
       invertButton.setBackground(Color.white);
       invertButton.setForeground(Color.black);
       _imgV.setInverted(false);
    }
    _inverted = !_inverted;
    _bigPanel.revalidate();
    _bigPanel.repaint();
    
  }

//==========================================
// get and set methods
//==========================================

   /**
    *   Returns the current value of the <em>dist</em> control.
    *   @return the current value of the <em>dist</em> control.
    */
   public double getDist() {
      double ret = 0;
      if(distSetter != null) {
         ret = ((Double)distSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   /**
    *   Sets the <em>dist</em> control with the given value.
    *   @param val the new value.
    */
   public void setDist(double val) {
      if(distSetter != null) {
         distSetter.setValue(val);
      }
   }
//......................................
   /**
    *   Returns the current value of the <em>pitch</em> control.
    *   @return the current value of the <em>pitch</em> control.
    */
   public double getPitch() {
      double ret = 0;
      if(pitchSetter != null) {
         ret = ((Double)pitchSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   /**
    *   Sets the <em>pitch</em> control with the given value.
    *   @param val the new value.
    */
   public void setPitch(double val) {
      if(pitchSetter != null) {
         pitchSetter.setValue(val);
      }
   }
//......................................
   /**
    *   Returns the current value of the <em>yaw</em> control.
    *   @return the current value of the <em>yaw</em> control.
    */
   public double getYaw() {
      double ret = 0;
      if(yawSetter != null) {
         ret = ((Double)yawSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   /**
    *   Sets the <em>yaw</em> control with the given value.
    *   @param val the new value.
    */
   public void setYaw(double val) {
      if(yawSetter != null) {
         yawSetter.setValue(val);
      }
   }
//......................................
   /**
    *   Returns the current value of the <em>roll</em> control.
    *   @return the current value of the <em>roll</em> control.
    */
   public double getRoll() {
      double ret = 0;
      if(rollSetter != null) {
         ret = ((Double)rollSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   /**
    *   Sets the <em>roll</em> control with the given value.
    *   @param val the new value.
    */
   public void setRoll(double val) {
      /* roll should not be adjustable if view mode is 'up_is_up' */
      if((rollSetter != null)&&(controlMenu_2_2.isSelected())) {
         rollSetter.setValue(val);
      }
   }
//......................................
   /**
    *   Returns the current value of the <em>zoom</em> control.
    *   @return the current value of the <em>zoom</em> control.
    */
   public int getZoom() {
      int ret = 0;
      if(zoomSetter != null) {
         ret = zoomSetter.getValue();
      }
      return ret;
   }

   /**
    *   Sets the <em>zoom</em> control with the given value.
    *   @param val the new value.
    */
   public void setZoom(int val) {
      if(zoomSetter != null) {
         zoomSetter.setValue(val);
      }
   }
//......................................
   /**
    *   Returns the current contents of the
    *   <em>mouse-click anatomy</em> text field.
    *   @return _anatomyStr.
    */
   public String getAnatomyStr() {
      return _anatomyStr;
   }

//......................................
   /**
    *   Returns the top level container of the SectionViewer.
    *   @return _contentPane.
    */
   public JPanel getContentPane() {
      return _contentPane;
   }

//......................................
   /**
    *   Returns the current contents of the
    *   <em>position feedback</em> text field.
    *   @return xyzTextField.
    */
   public JTextField getXyzTextField() {
      return xyzTextField;
   }

//......................................
   /**
    *   Returns the SectionViewer's JScrollPane.
    *   @return _imageScrollPane.
    */
   public JScrollPane getImageScrollPane() {
      return _imageScrollPane;
   }

//......................................
   /**
    *   Returns the object that contains the SectionViewer.
    *   <br>This will be either a JFrame or a JInternalFrame.
    *   @return _frame.
    */
  public Container getFrame() {
    return _frame;
  }

//......................................
   /**
    *   Sets the object that contains the SectionViewer.
    *   <br>This will be either a JFrame or a JInternalFrame.
    *   @param frame the object that contains the SectionViewer.
    */
  public void setFrame(Container frame) {
    _frame = frame;
  }

//......................................


//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
  /**
   *   Event handler for menu selections from 
   *   the <em>file</em> menu.
   */
  public class fileMenuHandler
      implements ActionListener {

    String currentDirStr;
    File currentDirFile;
    File selectedFile;
    JFileChooser chooser;

    public fileMenuHandler() {
      currentDirStr = System.getProperty("user.home");
      currentDirFile = new File(currentDirStr);
    }

    public void actionPerformed(ActionEvent e) {
      chooser = new JFileChooser(currentDirFile);
      if (e.getActionCommand() == fileMenu_1str) {
        // save image
	SVFileFilter filter = new SVFileFilter();
	chooser.setFileFilter(filter);
        int returnVal = chooser.showSaveDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          selectedFile = chooser.getSelectedFile();
          currentDirFile = chooser.getCurrentDirectory();
          saveImage(checkFilename(selectedFile.toString(), "jpg"));
        }
      }
      else if (e.getActionCommand() == fileMenu_2str) {
        // save view settings
        int returnVal = chooser.showSaveDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          selectedFile = chooser.getSelectedFile();
          currentDirFile = chooser.getCurrentDirectory();
          saveViewAsXML(checkFilename(selectedFile.toString(), "xml"));
        }
      }
      else if (e.getActionCommand() == fileMenu_3str) {
        // load view settings
        int returnVal = chooser.showOpenDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          selectedFile = chooser.getSelectedFile();
          currentDirFile = chooser.getCurrentDirectory();
          loadViewFromXML(selectedFile);
        }
      }
      else if (e.getActionCommand() == fileMenu_4str) {
        close();
      }
    } // actionPerformed()
  } // class fileMenuHandler

//---------------------------------------
  /**
   *   Event handler for menu selections from 
   *   the <em>control</em> menu.
   */
  public class controlMenuHandler
      implements ActionListener {

    double xa[] = new double[1];
    double ya[] = new double[1];
    double za[] = new double[1];

    double vals[] = null;

    PointEntry pe = null;
    JCheckBoxMenuItem src;

    public controlMenuHandler() {
    }

    public void actionPerformed(ActionEvent e) {

      // rotation ------------------------------
      if (e.getActionCommand().equals(controlMenu_1_1str)) {
        // see yaw pitch roll controls
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != controlMenu_1_1state) {
          if (src.isSelected()) {
            addStandardControls();
            _showStdCntrls = true;
          }
          else {
            removeStandardControls();
            _showStdCntrls = false;
            controlMenu_1_3state = false;
            controlMenu_1_3.setSelected(false);
          }
          controlMenu_1_1state = src.isSelected();
        }
      }
      else if (e.getActionCommand().equals(controlMenu_1_2str)) {
        // see user defined rotation controls
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != controlMenu_1_2state) {
          if (src.isSelected()) {
            addUserControls();
            _showUsrCntrls = true;
          }
          else {
            removeUserControls();
            _showUsrCntrls = false;
            controlMenu_1_3state = false;
            controlMenu_1_3.setSelected(false);
          }
          controlMenu_1_2state = src.isSelected();
        }
      }
      else if (e.getActionCommand().equals(controlMenu_1_3str)) {
        // see all rotation controls
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != controlMenu_1_3state) {
          if (src.isSelected()) {
	    addStandardControls();
	    _showStdCntrls = true;
            controlMenu_1_1state = true;
            controlMenu_1_1.setSelected(true);
	    addUserControls();
	    _showUsrCntrls = true;
            controlMenu_1_2state = true;
            controlMenu_1_2.setSelected(true);
          }
          else {
            removeStandardControls();
            _showStdCntrls = false;
            controlMenu_1_1state = false;
            controlMenu_1_1.setSelected(false);
            removeUserControls();
            _showUsrCntrls = false;
            controlMenu_1_2state = false;
            controlMenu_1_2.setSelected(false);
          }
          controlMenu_1_3state = src.isSelected();
        }
      } // view mode ------------------------------
      else if (e.getActionCommand().equals(controlMenu_2_1str)) {
        // up is up mode => roll is fixed
        setViewMode(controlMenu_2_1str);
      }
      else if (e.getActionCommand().equals(controlMenu_2_2str)) {
        setViewMode(controlMenu_2_2str);
      }
      else if (e.getActionCommand().equals(controlMenu_2_3str)) {
        /* this should be disabled normally
	   fixed line mode is set by the program */
        setViewMode(controlMenu_2_3str);
      }
      // fixed point ------------------------------
      else if (e.getActionCommand().equals(controlMenu_3_1str)) {
	 // change using mouse
	 _setAxisPoint = false;
	 _setFixedPoint = true;
	 setCursor(xhairCursor);
	 Runnable fpClick = new Runnable() {
	    public void run() {
	       if (!showMenu_4state) {
		  showMenu_4.doClick();
	       }
	    }
	 };
	 SwingUtilities.invokeLater(fpClick);
      }
      else if (e.getActionCommand().equals(controlMenu_3_2str)) {
        // change by typing coordinates
	 Runnable fpClick = new Runnable() {
	    public void run() {
	       if (!showMenu_4state) {
		  showMenu_4.doClick();
	       }
	    }
	 };
	 SwingUtilities.invokeLater(fpClick);
        _VSModel.getFixedPoint(xa, ya, za);
        vals = new double[3];
        vals[0] = xa[0];
        vals[1] = ya[0];
        vals[2] = za[0];
        _FPBounce = 0;
        pe = new PointEntry("Fixed Point Coordinates");
        pe.setSize(250, 250);
        pe.pack();
        pe.setVisible(true);
        pe.setInitVals(vals);
        pe.addChangeListener(new FixedPointAdaptor(pe));
        vals = null;
      }
      else if (e.getActionCommand().equals(controlMenu_3_3str)) {
        // reset
        resetFixedPoint(true);
      }
      // fixed line ------------------------------
      else if (e.getActionCommand().equals(controlMenu_4_1str)) {
	 // change using mouse
         if(!checkDist()) return; // check that dist = 0
	 _setFixedPoint = false;
	 _setAxisPoint = true;
	 setCursor(xhairCursor);
	 Runnable apClick = new Runnable() {
	    public void run() {
	       if (!showMenu_5state) {
		  showMenu_5.doClick();
	       }
	    }
	 };
	 SwingUtilities.invokeLater(apClick);
      }
      else if (e.getActionCommand().equals(controlMenu_4_2str)) {
        // change by typing coordinates
	boolean fpShowing = _fixedPoint;
	_fixedPoint = true;
	doShowFixedPoint();
	_fixedPoint = fpShowing;
	/*
	Runnable fpClick = new Runnable() {
	   public void run() {
	      if (!showMenu_4state) {
	         showMenu_4.doClick();
	      }
	   }
	};
	SwingUtilities.invokeLater(fpClick);
	*/
	Runnable apClick = new Runnable() {
	   public void run() {
	      if (!showMenu_5state) {
	         showMenu_5.doClick();
	      }
	   }
	};
	SwingUtilities.invokeLater(apClick);
        vals = _VSModel.getAxisPoint();
        _APBounce = 0;
        pe = new PointEntry("Fixed Line Coordinates");
        pe.setSize(250, 250);
        pe.pack();
        pe.setVisible(true);
        pe.setInitVals(vals);
        pe.addChangeListener(new FixedLineAdaptor(pe));
        vals = null;
      }
      else if (e.getActionCommand().equals(controlMenu_4_3str)) {
        // reset fixed line
        resetFixedLine(false);
        _VSModel.fireChange();
      }
      // reset ------------------------------
      else if (e.getActionCommand() == controlMenu_5str) {
        // reset user interface controls
        resetGUI();
        _VSModel.fireChange();
      }
    } // actionPerformed()

  } // class controlMenuHandler

//---------------------------------------
  /**
   *   Event handler for menu selections from 
   *   the <em>show</em> menu.
   */
  public class showMenuHandler
      implements ActionListener {

    JCheckBoxMenuItem src;

    public showMenuHandler() {
    }

    public void actionPerformed(ActionEvent e) {

      // cursor feedback ------------------------------
      if (e.getActionCommand().equals(showMenu_1str)) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != showMenu_1state) {
          if (src.isSelected()) {
            addFeedback();
          }
          else {
            removeFeedback();
          }
          showMenu_1state = src.isSelected();
        }
      }
      // intersection of views ------------------------------
      else if (e.getActionCommand().equals(showMenu_2str)) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != showMenu_2state) {
          if (src.isSelected()) {
            _showIntersection = true;
            doIntersection();
          }
          else {
            _showIntersection = false;
            removeIntersection();
          }
          showMenu_2state = src.isSelected();
        }
      }
      // mouse-click anatomy ------------------------------
      else if (e.getActionCommand().equals(showMenu_3str)) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != showMenu_3state) {
          if (src.isSelected()) {
            enableAnatomy();
            /* remove & disable thresholding */
            if (thresholdMenu_1state) {
              thresholdMenu_1.doClick();
            }
            if (thresholdMenu_2state) {
              thresholdMenu_2.doClick();
            }
            if (thresholdMenu_3state) {
              thresholdMenu_3.doClick();
            }
            thresholdMenu_4.doClick();
          }
          else {
            disableAnatomy();
            removeOverlay();
            resetAnatomyText();
          }
          showMenu_3state = src.isSelected();
        }
      }
      // fixed point ------------------------------
      else if (e.getActionCommand().equals(showMenu_4str)) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != showMenu_4state) {
          if (src.isSelected()) {
            _fixedPoint = true;
            doShowFixedPoint();
            _imgV.repaint();
          }
          else {
            _fixedPoint = false;
            removeFixedPoint();
          }
          showMenu_4state = src.isSelected();
        }
      }
      // fixed line ------------------------------
      else if (e.getActionCommand().equals(showMenu_5str)) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != showMenu_5state) {
          if (src.isSelected()) {
            _axisPoint = true;
            doShowFixedLine();
      revalidate();
      _bigPanel.repaint();
            //_imgV.repaint();
          }
          else {
            _axisPoint = false;
            removeAxisPoint();
          }
          showMenu_5state = src.isSelected();
        }
      }
    }
  } // class showMenuHandler
//---------------------------------------
  /**
   *   Event handler for menu selections from 
   *   the <em>help</em> menu.
   */
  public class helpMenuHandler
      implements ActionListener {

    // help system belongs to SectionViewer Utils class

    HelpBroker hb = null;
    JMenuItem mi = null;
    String IDStr = "";
    String viewStr = "TOC";

    public helpMenuHandler() {
    }

    public void actionPerformed(ActionEvent e) {
      if (null != _parent) {
        try {
          Method M1 = null;
          M1 = _parent.getClass().getMethod("getSVHelpBroker", null);
          hb = (HelpBroker) M1.invoke(_parent, null);
        }
        catch (InvocationTargetException ie) {
          System.out.println(ie.getMessage());
        }
        catch (NoSuchMethodException ne) {
          System.out.println("getSVHelpBroker: no such method");
          System.out.println(ne.getMessage());
        }
        catch (IllegalAccessException ae) {
          System.out.println(ae.getMessage());
        }
      }

      mi = (JMenuItem) e.getSource();

      try {
        if (mi.equals(helpMenu_1)) {
          viewStr = "TOC";
        }
        else if (mi.equals(helpMenu_2)) {
          viewStr = "Index";
        }
        else if (mi.equals(helpMenu_3)) {
          viewStr = "Search";
        }

        hb.setCurrentView(viewStr);
        hb.setDisplayed(true);
      }
      catch (Exception ex) {
        System.out.println(ex.getMessage());
      }

    } // actionPerformed()

  } // class helpMenuHandler

//---------------------------------------
  /**
   *   Event handler for menu selections from 
   *   the <em>threshold</em> menu.
   *   <p>Not in use.
   */
  public class thresholdMenuHandler
      implements ActionListener {

    /* the options are not part of a button group as they
       are not mutually exclusive */
    JCheckBoxMenuItem src;
    String message;
    int reply;

    public thresholdMenuHandler() {
    }

    public void actionPerformed(ActionEvent e) {

      if (e.getActionCommand() == thresholdMenu_1str) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != thresholdMenu_1state) {
          if (src.isSelected() == true) {
            /* allow a 'polygonal' constraint to be drawn */
            resetFeedbackText();
            enableThreshConstraint(true);
            enableThresholding(false);
            if (thresholdMenu_2state) {
              thresholdMenu_2state = false;
              thresholdMenu_2.setSelected(false);
            }
            /* remove & disable mouse click anatomy */
            if (showMenu_3state) {
              showMenu_3.doClick();
            }
          }
          else {
            enableThreshConstraint(false);
          }
          thresholdMenu_1state = src.isSelected();
        }
      }

      if (e.getActionCommand() == thresholdMenu_2str) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != thresholdMenu_2state) {
          if (src.isSelected() == true) {
            /* until bug is fixed vvvvvvvv */
            message = new String(
                "thresholding has a bug which may crash the program. Do you wish to proceed?");
            reply = JOptionPane.showConfirmDialog(
                null,
                message,
                "thresholding bug",
                JOptionPane.YES_NO_OPTION);
            if (reply == JOptionPane.YES_OPTION) {
              /* threshold based on grey value of clicked point */
              removeOverlay();
              resetFeedbackText();
              enableThresholding(true);
              if (thresholdMenu_1state) {
                thresholdMenu_1state = false;
                thresholdMenu_1.setSelected(false);
              }
              /* remove & disable mouse click anatomy */
              if (showMenu_3state) {
                showMenu_3.doClick();
              }
            }
            else {
              enableThresholding(false);
            }
            /* until bug is fixed ^^^^^^^ */
          }
          else {
            removeThreshold();
            enableThresholding(false);
          }
          thresholdMenu_2state = src.isSelected();
        }
      }

      if (e.getActionCommand() == thresholdMenu_3str) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != thresholdMenu_3state) {
          if (src.isSelected() == true) {
            // show 'polygonal' constraint
            _imgV.enableThreshConstraint(true);
            _imgV.repaint();
          }
          else {
            _imgV.enableThreshConstraint(false);
            _imgV.repaint();
          }
          thresholdMenu_3state = src.isSelected();
        }
      }

      if (e.getActionCommand() == thresholdMenu_4str) {
        // remove 'polygonal' constraint
        removeThreshConstraint();
        enableThreshConstraint(false);
        if (thresholdMenu_3state) {
          thresholdMenu_3state = false;
          thresholdMenu_3.setSelected(false);
        }
      }

    } // actionPerformed()

  } // class thresholdMenuHandler

//---------------------------------------
  /**
   *   Event handler for
   *   the <em>colour chooser</em> button.
   */
  public class planeColChooser
      implements ActionListener {

    public planeColChooser() {}

    public void actionPerformed(ActionEvent e) {
      Color c = JColorChooser.showDialog(null, "choose view colour",
                                         secColorClt.getBackground());
      if (null != c) {
        secColorClt.setBackground(c);
        updateIntersections();
      }
    }
  }

//---------------------------------------
  /**
   *   Event handler for
   *   the <em>invert</em> button.
   */
  public class invertButtonHandler implements ActionListener {

    public invertButtonHandler() {
    }

    public void actionPerformed(ActionEvent e) {
       // reverses the grey-level image
       invertSection();
    }
  }

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//---------------------------------------
// Control ADAPTORS, wsetter to wlz object
//---------------------------------------
  // change 'distance' using WSetter control
  /**
   *   Event handler for the <em>dist</em> slider control.
   */
  public class WSetterToDistAdaptor
      implements ChangeListener {

    WSetter control;
    ViewStructModel VSmodel;
    WlzObjModel OBJmodel;
    Vector vec = null;
    double val;

    public WSetterToDistAdaptor(WSetter cntrl,
                                ViewStructModel mdl1,
                                WlzObjModel mdl2) {
      VSmodel = mdl1;
      OBJmodel = mdl2;
      control = cntrl;
    }

    public void stateChanged(ChangeEvent e) {
      vec = control.getValue();
      switch (control.getType()) {
        case INTEGER:
          Integer iVal = (Integer) vec.elementAt(0);
          val = (double) iVal.intValue();
          break;
        case FLOAT:
          Float fVal = (Float) vec.elementAt(0);
          val = (double) fVal.floatValue();
          break;
        case DOUBLE:
          Double dVal = (Double) vec.elementAt(0);
          val = dVal.doubleValue();
          break;
        default:
          break;
      }
      VSmodel.setDist(val);
      if((val != 0.0) &&
         !VSmodel.getViewMode().equals("WLZ_FIXED_LINE_MODE")) {
         setViewMode(controlMenu_2_1str); // up is up
      }
    }
  }

//---------------------------------------
  // change 'pitch' using WSetter control
  /**
   *   Event handler for the <em>pitch</em> slider control.
   */
  public class WSetterToPitchAdaptor
      implements ChangeListener {

    WSetter control;
    ViewStructModel VSmodel;
    WlzObjModel OBJmodel;
    double valArr[] = new double[1];
    Vector vec = null;
    double val;

    public WSetterToPitchAdaptor(WSetter cntrl, ViewStructModel mdl1,
                                 WlzObjModel mdl2) {
      VSmodel = mdl1;
      OBJmodel = mdl2;
      control = cntrl;
    }

    public void stateChanged(ChangeEvent e) {
       vec = control.getValue();
       switch (control.getType()) {
	  case INTEGER:
	     Integer iVal = (Integer) vec.elementAt(0);
	     val = (double) iVal.intValue();
	     break;
	  case FLOAT:
	     Float fVal = (Float) vec.elementAt(0);
	     val = (double) fVal.floatValue();
	     break;
	  case DOUBLE:
	     Double dVal = (Double) vec.elementAt(0);
	     val = dVal.doubleValue();
	     break;
	  default:
	     break;
       }
       VSmodel.setPhiDeg(val);
       VSmodel.getZeta(valArr);
       val = (180.0 / Math.PI) * valArr[0];
       rollSetter.setSliderEnabled(true);
       rollSetter.setValue(val);
       rollSetter.setSliderEnabled(false);
       if (!control.isAdjusting()) {
	  VSmodel.getDist(valArr);
	  setDistLimits(valArr[0]);
       }
    }
  }

//---------------------------------------
  // change 'yaw' using WSetter control
  /**
   *   Event handler for the <em>yaw</em> slider control.
   */
  public class WSetterToYawAdaptor
      implements ChangeListener {

    WSetter control;
    ViewStructModel VSmodel;
    WlzObjModel OBJmodel;
    double valArr[] = new double[1];
    Vector vec = null;
    double val;

    public WSetterToYawAdaptor(WSetter cntrl,
                               ViewStructModel mdl1,
                               WlzObjModel mdl2) {
      VSmodel = mdl1;
      OBJmodel = mdl2;
      control = cntrl;
    }

    public void stateChanged(ChangeEvent e) {
      vec = control.getValue();
      switch (control.getType()) {
        case INTEGER:
          Integer iVal = (Integer) vec.elementAt(0);
          val = (double) iVal.intValue();
          break;
        case FLOAT:
          Float fVal = (Float) vec.elementAt(0);
          val = (double) fVal.floatValue();
          break;
        case DOUBLE:
          Double dVal = (Double) vec.elementAt(0);
          val = dVal.doubleValue();
          break;
        default:
          break;
      }
      VSmodel.setThetaDeg(val);
      if (!control.isAdjusting()) {
        VSmodel.getDist(valArr);
        setDistLimits(valArr[0]);
      }
    }
  }

//---------------------------------------
  // change 'roll' using WSetter control
  /**
   *   Event handler for the <em>roll</em> slider control.
   */
  public class WSetterToRollAdaptor
      implements ChangeListener {

    WSetter control;
    ViewStructModel VSmodel;
    WlzObjModel OBJmodel;
    double valArr[] = new double[1];
    Vector vec = null;
    double val;

    public WSetterToRollAdaptor(WSetter cntrl,
                                ViewStructModel mdl1,
                                WlzObjModel mdl2) {
      VSmodel = mdl1;
      OBJmodel = mdl2;
      control = cntrl;
    }

    public void stateChanged(ChangeEvent e) {
      vec = control.getValue();
      switch (control.getType()) {
        case INTEGER:
          Integer iVal = (Integer) vec.elementAt(0);
          val = (double) iVal.intValue();
          break;
        case FLOAT:
          Float fVal = (Float) vec.elementAt(0);
          val = (double) fVal.floatValue();
          break;
        case DOUBLE:
          Double dVal = (Double) vec.elementAt(0);
          val = dVal.doubleValue();
          break;
        default:
          break;
      }
      VSmodel.setZetaDeg(val);
      if (!control.isAdjusting()) {
        VSmodel.getDist(valArr);
        setDistLimits(valArr[0]);
      }
    }
  }
//---------------------------------------
  // change fixed line rotation using WSetter control
  /**
   *   Event handler for the <em>fixed line</em> slider control.
   */
  public class WSetterToFixedRotAdaptor
     implements ChangeListener {

        WSetter control;
        ViewStructModel VSmodel;
        WlzObjModel OBJmodel;
        double valArr[] = new double[1];
        Vector vec = null;
        double val;

        public WSetterToFixedRotAdaptor(WSetter cntrl, ViewStructModel mdl1,
              WlzObjModel mdl2) {
           VSmodel = mdl1;
           OBJmodel = mdl2;
           control = cntrl;
        }

        public void stateChanged(ChangeEvent e) {
           if(!_fixedLineRotation) return;
           vec = control.getValue();
           switch (control.getType()) {
              case INTEGER:
                 Integer iVal = (Integer) vec.elementAt(0);
                 val = (double) iVal.intValue();
                 break;
              case FLOAT:
                 Float fVal = (Float) vec.elementAt(0);
                 val = (double) fVal.floatValue();
                 break;
              case DOUBLE:
                 Double dVal = (Double) vec.elementAt(0);
                 val = dVal.doubleValue();
                 break;
              default:
                 break;
           }

           VSmodel.setPsiDeg(val);

           /* adjust the Yaw Pitch and Roll sliders */
           //if (!control.isAdjusting()) {
              VSmodel.getTheta(valArr);
              double dangle = (180.0 / Math.PI) * valArr[0];
              yawSetter.setSliderEnabled(true);
              yawSetter.setValue(dangle);
              yawSetter.setSliderEnabled(false);

              VSmodel.getPhi(valArr);
              dangle = (180.0 / Math.PI) * valArr[0];
              pitchSetter.setSliderEnabled(true);
              pitchSetter.setValue(dangle);
              pitchSetter.setSliderEnabled(false);

              VSmodel.getZeta(valArr);
              dangle = (180.0 / Math.PI) * valArr[0];
              rollSetter.setSliderEnabled(true);
              rollSetter.setValue(dangle);
              rollSetter.setSliderEnabled(false);

              VSmodel.getDist(valArr);
              setDistLimits(valArr[0]);
           //}
        }
     }



//---------------------------------------
// View ADAPTORS
//---------------------------------------
  // change the displayed section when the view structure changes
  /**
   *   Event handler that changes the current section when
   *   any aspect of the <em>View Structure</em> changes.
   */
  public class ViewStructModelToViewAdaptor
      implements ChangeListener {

    ViewStructModel VSmodel;
    WlzThreeDViewStruct VS;
    WlzObjModel OBJmodel;
    WlzImgView view;
    Dimension imgSize;
    Dimension newSize;
    double mag;

    WlzObject section;
    WlzObject oSection;
    WlzObject oObj;
    WlzIBox2 bbox2;
    Vector anatVec;
    WlzObject obj2D[] = null;

    public ViewStructModelToViewAdaptor(ViewStructModel mdl1,
                                        WlzObjModel mdl2,
                                        WlzImgView vw) {
      VSmodel = mdl1;
      OBJmodel = mdl2;
      view = vw;
      VS = VSmodel.getViewStruct();
    }

    public void stateChanged(ChangeEvent e) {
      if (_debug) {
        System.out.println("enter ViewStructModelToViewAdaptor");
      }

      if (_enabled == false) return;

      section = null;
      oObj = null;
      oSection = null;

      int num = AnatKey.getNRows();

      obj2D = new WlzObject[num];
      boolean viz[] = new boolean[num];
      Color col[] = new Color[num];
      
      for (int i = 0; i < num; i++) {
        obj2D[i] = null;
        viz[i] = false;
	col[i] = null;
      }

      view.clearThreshold();
      view.enableThreshConstraint(false);

      resetPosGreyText();
      section = OBJmodel.getSection();

//-------------------------
      // draw the grey level
      try {
        view.setWlzObj(section);
      }
      catch (WlzException e1) {
        System.out.println("ViewStructModelToViewAdaptor #1");
        System.out.println(e1.getMessage());
      }
//-------------------------
      doShowFixedPoint();
      doShowFixedLine();
//-------------------------
      // draw the anatomy components
      try {
        if (null != _parent) {
          try {
            Method M1 = null;
            M1 = _parent.getClass().getMethod("getAnatomyElements", null);
            anatVec = (Vector)M1.invoke(_parent, null);
          }
          catch (InvocationTargetException ie) {
            System.out.println(ie.getMessage());
          }
          catch (NoSuchMethodException ne) {
            System.out.println("getAnatomyElements: no such method");
            System.out.println(ne.getMessage());
          }
          catch (IllegalAccessException ae) {
            System.out.println(ae.getMessage());
          }
        }

        if ( (anatVec != null) && (viz != null) && (obj2D != null)) {
	  Enumeration els = anatVec.elements();

	  int i = 0;
	  AnatomyElement el = null;

	  while(els.hasMoreElements()) {
	    el = (AnatomyElement)els.nextElement();
	    if (el != null) {
	      viz[i] = el.isVisible();
	      col[i] = AnatKey.getColor(el.getIndx());
	      obj2D[i] = _OBJModel.makeSection(
		  el.getObj(),
		  VS);
	      if (obj2D[i] != null) {
		if (WlzObject.WlzBndObjGetType(obj2D[i]) ==
		    WlzObjectType.WLZ_EMPTY_OBJ) {
		  obj2D[i] = null;
		}
	      }
	    }
	    else {
	      obj2D[i] = null;
	      viz[i] = false;
	      col[i] = null;
	    }
	    i++;
	  } // while
          view.setAnatomyObj(obj2D, viz, col);
        }
      }
      catch (WlzException e2) {
        System.out.println("ViewStructModelToViewAdaptor #2");
        System.out.println(e2.getMessage());
      }
//-------------------------
      // draw the 'mouse-click' anatomy
      if (_thresholding != true) {
        view.clearOverlay();
        if ( (_maybeObj != null) && (_maybeObj.empty() == false)) {
          oObj = (WlzObject) _maybeObj.peek();
          if (oObj != null) {
            oSection = OBJmodel.makeSection(oObj, VS);
          }
          if (oSection != null) {
            // trap WLZ_ERR_OBJECT_TYPE to see
            // if section 'misses' anatomy component
            try {
              view.setOverlayObj(oSection);
            }
            catch (WlzException e3) {
              System.out.println("ViewStructModelToViewAdaptor #3");
              System.out.println(e3.getMessage());
              if ( (e3.getMessage()).equals("WLZ_ERR_OBJECT_TYPE")) {
                System.out.println("oSection misses anatomy");
              }
            }
          }
          else {
            //System.out.println("oSection == null");
          }
        } // if((_maybeObj != null) && ...
      } // if(_thresholding != true)
//-------------------------
      imgSize = view.getImgSize();
      mag = view.getMag();
      newSize = new Dimension(
          (int) (imgSize.width * mag + 10),
          (int) (imgSize.height * mag + 10));
      setViewTitle();
      _bigPanel.setPreferredSize(newSize);
      revalidate();
      _bigPanel.repaint();

      updateIntersections();

      if (_debug) {
        System.out.println("exit ViewStructModelToViewAdaptor");
      }

    } // stateChanged

  } // ViewStructModelToViewAdaptor

//---------------------------------------
  // connect mouse events to display of x-y position
  /**
   *   Event handler that updates the 
   *   <em>position</em> feedback text field when
   *   _imgV fires an event. 
   */
  public class WlzImgViewToPosAdaptor
      implements ChangeListener {

    WlzImgView ImgModel;
    WlzObjModel ObjModel;
    ViewStructModel VSModel;
    WlzThreeDViewStruct VS;
    JTextField view;
    Point pos;
    double pt3d[];

    public WlzImgViewToPosAdaptor(WlzImgView mdl1,
                                  WlzObjModel mdl2,
                                  ViewStructModel mdl3,
                                  JTextField tfpos) {
      ImgModel = mdl1;
      ObjModel = mdl2;
      VSModel = mdl3;
      view = tfpos;
      VS = VSModel.getViewStruct();
    }

    public void stateChanged(ChangeEvent e) {
      pos = ImgModel.getPos();
      pt3d = ObjModel.get3DPoint(pos, VS);
      view.setText(Math.round(pt3d[0]) + ", " +
                   Math.round(pt3d[1]) + ", " +
                   Math.round(pt3d[2]));
    }

  }

//---------------------------------------
  // connects mouse events to display of grey level value
  /**
   *   Event handler that updates the 
   *   <em>grey value</em> feedback text field when
   *   _imgV fires an event. 
   */
  public class WlzImgViewToGreyAdaptor
      implements ChangeListener {

    WlzImgView model;
    JTextField view;
    int val;

    public WlzImgViewToGreyAdaptor(WlzImgView mdl, JTextField pos) {
      model = mdl;
      view = pos;
    }

    public void stateChanged(ChangeEvent e) {
      val = model.getGreyVal();
      view.setText(Integer.toString(val));
    }

  }

//---------------------------------------
  // connects mouse events to display of anatomy
  /**
   *   Event handler that updates the 
   *   <em>mouse-click anatomy</em> feedback text field when
   *   _imgV fires an event. 
   */
  public class WlzImgViewToAnatAdaptor
      implements ChangeListener {

    ViewStructModel VSmodel;
    WlzThreeDViewStruct VS;
    WlzObjModel OBJmodel;
    WlzImgView IMGmodel;
    JTextField view;
    Point xypos;

    WlzObject oSection = null;
    WlzObject oObj = null;

    double pt3d[];

    public WlzImgViewToAnatAdaptor(ViewStructModel mdl1,
                                   WlzObjModel mdl2,
                                   WlzImgView mdl3,
                                   ScrollableTextField anat) {
//                                   JTextField anat) {
      VSmodel = mdl1;
      VS = VSmodel.getViewStruct();
      OBJmodel = mdl2;
      IMGmodel = mdl3;
      view = anat;
    }

    public void stateChanged(ChangeEvent e) {

      oObj = null;
      oSection = null;
      if ( (_thresholding == false) &&
          (_threshConstraint == false) &&
          (_anatomy == true)) {
        xypos = IMGmodel.getPos();
        pt3d = OBJmodel.get3DPoint(new Point(xypos.x, xypos.y), VS);
        /*
           System.out.println("2d pt: "+xypos.x+", "+xypos.y);
           System.out.println("3d pt: "+pt3d[0]+", "+pt3d[1]+", "+pt3d[2]);
         */
        _anatomyStr = getAnatomyAt(pt3d);
        view.setText(_anatomyStr);

        // overlay component onto greyscale image
        if ( (_maybeObj != null) && (_maybeObj.empty() == false)) {
          oObj = (WlzObject) _maybeObj.peek();
          //printFacts(oObj);
          if (oObj != null) {
            oSection = OBJmodel.makeSection(oObj, VS);
          }
          else {
            IMGmodel.clearOverlay();
          }
          if (oSection != null) {
            //printFacts(oSection);
            try {
              IMGmodel.setOverlayObj(oSection);
            }
            catch (WlzException err) {
              System.out.println("WlzException #5");
              System.out.println(err.getMessage());
            }
          }
          else {
            IMGmodel.clearOverlay();
          }
        }
        else {
          view.setText("");
          IMGmodel.clearOverlay();
        }
      }
      else {
        view.setText("");
        IMGmodel.clearOverlay();
      }
      IMGmodel.repaint();
    }
  }

//---------------------------------------
  // change 'zoom' using Zoom control
  /**
   *   Event handler for
   *   the <em>zoom</em> control.
   */
  public class ZoomToZoomAdaptor
      implements ChangeListener {

    Zoom control;
    WlzImgView view;
    double zfactor;
    double mag;
    int zval;
    int min, max, val, ext;
    int imgW, imgH;

    Dimension imgSize;
    Dimension newSize;

    BoundedRangeModel BRM = null;

    public ZoomToZoomAdaptor(Zoom cntrl,
                                WlzImgView vw) {
      control = cntrl;
      view = vw;
    }

    public void stateChanged(ChangeEvent e) {
      //System.out.println("ZoomToZoomAdaptor: stateChanged");
      if (view != null) {
        zval = control.getValue();
        zfactor = (double)zval;
        mag = zf2mag(zfactor);
        imgSize = view.getImgSize();
        view.setMag(mag);
        imgW = (int) (imgSize.width * mag);
        imgH = (int) (imgSize.height * mag);
        newSize = new Dimension(imgW+10, imgH+10);

        _bigPanel.setPreferredSize(newSize);
        _bigPanel.revalidate();
        _bigPanel.repaint();

      }
    }
  }

//---------------------------------------
  // connects mouse events to feedback display
  /**
   *   Event handler that updates <em>_imgV</em> with
   *   the current cursor position in the 2D grey-level image,
   *   then causes <em>_imgV</em> to fire a changeEvent.
   */
  public class MouseToFBAdaptor
      implements MouseListener, MouseMotionListener {

    WlzImgView model;
    Point xypos;
    double mag;

    Point pt;

    public MouseToFBAdaptor(WlzImgView mdl) {
      model = mdl;
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mousePressed(MouseEvent e) {
      mag = model.getMag();
      double xnorm = (e.getX()) / mag;
      double ynorm = (e.getY()) / mag;
      model.updateStats(xnorm, ynorm);
      model.fireChange();
    }

    public void mouseReleased(MouseEvent e) {}

    public void mouseClicked(MouseEvent e) {
    }

    public void mouseMoved(MouseEvent e) {
      if (e.getModifiers() == e.SHIFT_MASK) {
        mag = model.getMag();
        double xnorm = (e.getX()) / mag;
        double ynorm = (e.getY()) / mag;
        model.updateStats(xnorm, ynorm);
        model.fireChange();
      }
    }

    public void mouseDragged(MouseEvent e) {
      mag = model.getMag();
      double xnorm = (e.getX()) / mag;
      double ynorm = (e.getY()) / mag;
      model.updateStats(xnorm, ynorm);
      model.fireChange();
    }
  }

//---------------------------------------
  // connects mouse events to thresholding
  /**
   *   Event handler that <em>thresholds</em> a region of
   *   the 2D grey-level image centred on an initial mouse press
   *   and controlled by the amount of mouse drag left or right.
   */
  public class MouseToThresholdAdaptor
      implements MouseListener, MouseMotionListener {

    WlzImgView view;
    WlzObjModel model;
    WlzObject section = null;
    WlzObject section1 = null;
    WlzObject section2 = null;
    WlzObject section3 = null;
    WlzObject section4 = null;
    WlzObject selectedObj = null;
    WlzDBox2 bbox2 = null;

    Dimension imgSize;
    Dimension newSize;

    double mag;
    int val;
    int inc = 0;
    int newval = 0;
    int hi = WlzThresholdType.WLZ_THRESH_HIGH;
    int lo = WlzThresholdType.WLZ_THRESH_LOW;

    double xorg;
    double yorg;
    double xpt;

    //Point pt2d = new Point();
    Point.Double pt2d = new Point.Double();

    public MouseToThresholdAdaptor(WlzImgView vw,
                                   WlzObjModel mdl) {
      view = vw;
      model = mdl;
    }

    public void mousePressed(MouseEvent e) {
      if ( (!_anatomy) && (_thresholding)) {
        view.clearThreshold();
        view.repaint();
        section = model.getSection();
        try {
          bbox2 = section.WlzBoundingBox2D(section);
        }
        catch (WlzException we) {
          System.out.println("WlzException #?");
          System.out.println(we.getMessage());
        }
        section1 = model.smooth(section, 1.0, 1.0, 0, 0);
        if (section1 != null) {
          mag = view.getMag();
          xorg = (e.getX()) / mag;
          yorg = (e.getY()) / mag;
          view.updateStats(xorg, yorg);
          val = view.getGreyVal();
          if (bbox2 != null) {
            pt2d.setLocation(bbox2.xMin + xorg, bbox2.yMin + yorg);
          }
          if (_threshConstraint) {
            section2 = model.constrain(section1, _constraintObj);
          }
          else {
            section2 = section1;
          }
        }
      }
    }

    public void mouseDragged(MouseEvent e) {
      /*
              if mouse dragged to right increase threshold window
              if mouse dragged to left decrease threshold window to min = 0
       */
      selectedObj = null;
      if ( (!_anatomy) && (_thresholding)) {
        if (section2 != null) {
          mag = view.getMag();
          xpt = (e.getX()) / mag;
          /*
                System.out.print("mouse at point: ");
                System.out.println((e.getX())/mag+", "+(e.getY())/mag);
                System.out.println("change = "+ (int)change(xpt));
           */
          int ch = (int) change(xpt);
          inc = ch > 0 ? ch : 0;
          if (inc > 0) {
            /* inc here is position not grey value !!!!!
               it is just a handy way to increase threshold */
            newval = (val + inc) > 255 ? 255 : (val + inc);
            section3 = model.threshold(section2, newval, lo);
            if (section3 != null) {
              newval = (val - inc) > 0 ? (val - inc) : 0;
              section4 = model.threshold(section3, newval, hi);
            }
            if (section4 != null) {
              // call native method to select domain containing point
              selectedObj = model.select(section4,
                                         pt2d.getX(),
                                         pt2d.getY());
              if (selectedObj != null) {
                try {
                  view.setThresholdObj(selectedObj);
                }
                catch (WlzException err) {
                  System.out.println("WlzException #7");
                  System.out.println(err.getMessage());
                }
              }
            }
            view.repaint();
          }
          else {
            view.clearThreshold();
          }
        }
      }
    }

    public void mouseReleased(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseClicked(MouseEvent e) {
    }

    public void mouseMoved(MouseEvent e) {
    }

    protected double change(double x) {
      return (x - xorg);
    }
  }

//---------------------------------------
  // connects mouse events to threshConstraint
  /**
   *   Event handler that draws a <em>threshold constraint</em>
   *   contour, starting at an initial mouse press, finishing when
   *   the mouse is released and defined by the intervening mouse drag.
   */
  public class MouseToThreshConstraintAdaptor
      implements MouseListener, MouseMotionListener {

    WlzImgView view;
    WlzObjModel model;
    WlzObject section = null;
    WlzIBox2 bbox2 = null;

    Vector xpts;
    Vector ypts;

    int npts = 0;

    double mag;
    double xorg;
    double yorg;
    double x;
    double y;

    Point pt2d = new Point();

    public MouseToThreshConstraintAdaptor(WlzImgView vw,
                                          WlzObjModel mdl) {
      view = vw;
      model = mdl;
    }

    public void mousePressed(MouseEvent e) {
      if ( (_anatomy == false) &&
          (_thresholding == false) &&
          (_threshConstraint == true)) {
        xpts = new Vector(100, 10);
        ypts = new Vector(100, 10);
        section = model.getSection();
        if (section != null) {
          try {
            bbox2 = section.WlzBoundingBox2I(section);
          }
          catch (WlzException we) {
            System.out.println("WlzException #8");
            System.out.println(we.getMessage());
          }
          mag = view.getMag();
          xorg = (e.getX()) / mag;
          yorg = (e.getY()) / mag;
          view.updateStats(xorg, yorg);
          view.clearThreshConstraint();
          view.enableThreshConstraint(true);
          xpts.add(new Float(xorg));
          ypts.add(new Float(yorg));
          view.repaint();
          /*
                 System.out.println("mouse pressed at");
                 System.out.print("xorg,yorg "+xorg+", "+yorg+" (");
                 System.out.println(pt2d.getX()+", "+pt2d.getY()+")");
           */
        }
      }
    }

    public void mouseDragged(MouseEvent e) {
      if ( (_anatomy == false) &&
          (_thresholding == false) &&
          (_threshConstraint == true)) {
        x = (e.getX()) / mag;
        y = (e.getY()) / mag;
        xpts.add(new Float(x));
        ypts.add(new Float(y));
        view.setThreshConstraint(xpts, ypts);
        view.repaint();
      }
    }

    public void mouseReleased(MouseEvent e) {
      if ( (_anatomy == false) &&
          (_thresholding == false) &&
          (_threshConstraint == true)) {

        view.closeThreshConstraint();
        view.repaint();
        constraintToWlz(bbox2);

      }
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseClicked(MouseEvent e) {
    }

    public void mouseMoved(MouseEvent e) {
    }

  }

//---------------------------------------
  // connects mouse entered events to focus notification
  /**
   *   Event handler which fires an ActionEvent
   *   when the mouse is clicked inside a SectionViewer.
   */
  public class FocusAdaptor
      implements MouseListener {

    ActionEvent event = null;
    String actionCommand = "focus";

    public FocusAdaptor() {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseClicked(MouseEvent e) {
      event = new ActionEvent(e.getSource(),
                              ActionEvent.ACTION_PERFORMED,
                              actionCommand);
      fireEvent(event);
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mousePressed(MouseEvent e) {
    }

    public void mouseReleased(MouseEvent e) {
    }

  }

//---------------------------------------
  // connects mouse events to fixed point / fixed line setting
  /**
   *   Event handler which sets the <em>Fixed Point</em> or 
   *   <em>2nd Fixed Point</em> when the mouse button is released.
   */
  public class MouseToFPAdaptor
      implements MouseListener, MouseMotionListener {

    double mag;
    WlzImgView ImgModel;
    WlzObjModel ObjModel;
    ViewStructModel VSModel;
    WlzThreeDViewStruct VS;
    Point pos;
    double FP3D[] = null;
    double FP2D[] = null;
    boolean fpShowing = false;

    public MouseToFPAdaptor(WlzImgView mdl1,
                            WlzObjModel mdl2,
                            ViewStructModel mdl3) {

      ImgModel = mdl1;
      ObjModel = mdl2;
      VSModel = mdl3;
      VS = VSModel.getViewStruct();
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
      if (!_setFixedPoint && !_setAxisPoint) return;
      fpShowing = _fixedPoint;
      _fixedPoint = true;
      doShowFixedPoint();
      _fixedPoint = fpShowing;
    }

    public void mousePressed(MouseEvent e) {
      if (!_setFixedPoint && !_setAxisPoint) return;
      pos = ImgModel.getPos();
      FP3D = ObjModel.get3DPoint(pos, VS);
      // restrict point to inside image
      clipPoint3D(FP3D);
      FP2D = ObjModel.get2DPoint(FP3D, VS);
      if(_setFixedPoint) {
	 _imgV.setFixedPointVec(FP2D);
      } else if(_setAxisPoint) {
	 _imgV.setAxisPointArr(FP2D);
      }
    }

    public void mouseReleased(MouseEvent e) {

      if (!_setFixedPoint && !_setAxisPoint) return;
      pos = ImgModel.getPos();
      FP3D = ObjModel.get3DPoint(pos, VS);
      // restrict point to inside image
      clipPoint3D(FP3D);

      if(_setFixedPoint) {
	 VSModel.setFixedPoint(FP3D[0],
			       FP3D[1],
			       FP3D[2]);
	 _setFixedPoint = false;
      } else if(_setAxisPoint) {
	 VSModel.setAxisPoint(FP3D[0],
	                      FP3D[1],
			      FP3D[2]);
	 setFixedLineParams();
	 setViewMode("fixed line");
	 _setAxisPoint = false;
	 if(!fpShowing) {
	    removeFixedPoint();
	 }
      }

      setDistLimits(0.0);
      setCursor(defCursor);
    }

    public void mouseClicked(MouseEvent e) {
    }

    public void mouseMoved(MouseEvent e) {
    }

    public void mouseDragged(MouseEvent e) {
      if (!_setFixedPoint && !_setAxisPoint) return;
      pos = ImgModel.getPos();
      FP3D = ObjModel.get3DPoint(pos, VS);
      // restrict point to inside image bounding box
      clipPoint3D(FP3D);
      FP2D = ObjModel.get2DPoint(FP3D, VS);

      if(_setFixedPoint) {
	 _imgV.setFixedPointVec(FP2D);
      } else if(_setAxisPoint) {
	 _imgV.setAxisPointArr(FP2D);
      }

    }
  }
//---------------------------------------
  /**
   *   Event handler which sets the <em>Fixed Point</em> 
   *   to values obtained from a <em>PointEntry</em> dialogue.
   */
  public class FixedPointAdaptor
      implements ChangeListener {

    double vals[];
    PointEntry _pe;

    public FixedPointAdaptor(PointEntry pe) {
      _pe = pe;
    }

    public void stateChanged(ChangeEvent e) {
      if (_FPBounce > 0) {
        vals = _pe.getValues();
        _VSModel.setFixedPoint(vals[0],
                               vals[1],
                               vals[2]);
        setDistLimits(0.0);
        Runnable fpClick = new Runnable() {
          public void run() {
            if (!showMenu_4state)
              showMenu_4.doClick();
          }
        };
        SwingUtilities.invokeLater(fpClick);
      }
      else {
        _FPBounce++;
      }
    }
  }

//-------------------------------------------------------------
  /**
   *   Event handler which sets the <em>2nd Fixed Point</em> 
   *   to values obtained from a <em>PointEntry</em> dialogue.
   */
  public class FixedLineAdaptor
      implements ChangeListener {

    double vals[];
    PointEntry _pe;

    public FixedLineAdaptor(PointEntry pe) {
      _pe = pe;
    }

    public void stateChanged(ChangeEvent e) {
      if (_APBounce > 0) {
        vals = _pe.getValues();
        _VSModel.setAxisPoint(vals[0],
                             vals[1],
                             vals[2]);

        setFixedLineParams();
        setViewMode("fixed line");
        _setAxisPoint = false;
        if(!_fixedPoint) {
           removeFixedPoint();
        }
        setDistLimits(0.0);
      }
      else {
        _APBounce++;
      }
    }
  }

//=============================================================
//-------------------------------------------------------------
// thread stuff
//-------------------------------------------------------------

  //Object _lock = new Object();
//-------------------------------------------------------------
  // stack
  /**
   *   Thread which is responsible for copying stacks of Woolz objects 
   *   and their filenames from _anatBuilder.
   */
  Runnable rStack = new Runnable() {
    public void run() {
      //System.out.println("stack thread started ");
      synchronized (_anatBuilder._lock) {
        while (_anatBuilder.getDone() == false) {
          try {
            _anatBuilder._lock.wait();
          }
          catch (InterruptedException iex2) {
            System.out.println("Stack thread interrupted2 ");
            System.out.println(iex2.getMessage());
          }
        }

        try {
          //System.out.println("stack thread working");
          if (_anatBuilder != null) {
            _embFileStack = (Stack) _anatBuilder.getEmbFileStack().clone();
            _xembFileStack = (Stack) _anatBuilder.getXEmbFileStack().clone();
            _embObjStack = (Stack) _anatBuilder.getEmbObjStack().clone();
            _xembObjStack = (Stack) _anatBuilder.getXEmbObjStack().clone();
            showMenu_3.setEnabled(true);
          }
          //System.out.println("stack thread finished");
        }
        catch (Exception ex) {
          System.out.println("SectionViewer stack thread");
          ex.printStackTrace();
        }
      } // sync1
    }
  };

//=============================================================
//-------------------------------------------------------------
// ChangeEvents
//-------------------------------------------------------------

//-------------------------------------------------------------
// handle all _objects that are interested in changes
//-------------------------------------------------------------
  // keep track of all the listeners to this 'model'
  /**
   *   A list of ChangeListeners which are
   *   listening for events fired from the SectionViewer.
   */
  protected EventListenerList changeListeners =
      new EventListenerList();

//-------------------------------------------------------------
  // add a listener to the register
  /**
   *   Adds a ChangeListener to the EventListenerList.
   */
  public void addChangeListener(ChangeListener x) {
    changeListeners.add(ChangeListener.class, x);

    // bring it up to date with current state
    x.stateChanged(new ChangeEvent(this));
  }

//-------------------------------------------------------------
  // remove a listener from the register
  /**
   *   Removes a ChangeListener from the EventListenerList.
   */
  public void removeChangeListener(ChangeListener x) {
    changeListeners.remove(ChangeListener.class, x);
  }

//-------------------------------------------------------------
  /**   An event that will be fired from SectionViewer. */
  private ChangeEvent ce;

  /**  A local copy of the list of ChangeListeners. */
  private Object[] listeners;

  /**  One of the list of Change Listeners. */
  private ChangeListener cl;

  /**
   *   Fires a ChangeEvent from the SectionViewer.
   */
  protected void fireChange() {
    if (_enabled == true) {
      //System.out.println("firing an event from SectionViewer");
      // Create the event:
      ce = new ChangeEvent(this);
      // Get the listener list
      listeners = changeListeners.getListenerList();
      // Process the listeners last to first
      // List is in pairs, Class and instance
      for (int i
           = listeners.length - 2; i >= 0; i -= 2) {
        if (listeners[i] == ChangeListener.class) {
          cl = (ChangeListener) listeners[i + 1];
          cl.stateChanged(ce);
        }
      }
    }
  } // fireChange

//=============================================================
//-------------------------------------------------------------
// ActionEvents
//-------------------------------------------------------------

//-------------------------------------------------------------
  // keep track of all the listeners to this 'model'
  /**
   *   A list of ActionListeners which are
   *   listening for events fired from the SectionViewer.
   */
  protected EventListenerList actionListeners =
      new EventListenerList();

//-------------------------------------------------------------
  // add a listener to the register
  /**
   *   Adds an ActionListener to the EventListenerList.
   */
  public void addActionListener(ActionListener x) {
    actionListeners.add(ActionListener.class, x);
  }

//-------------------------------------------------------------
  // remove a listener from the register
  /**
   *   Removes an ActionListener from the EventListenerList.
   */
  public void removeActionListener(ActionListener x) {
    actionListeners.remove(ActionListener.class, x);
  }

//-------------------------------------------------------------
  /**  A local copy of the list of ActionListeners. */
  private Object[] listeners2;

  /**  One of the list of Action Listeners. */
  private ActionListener al;

  /**
   *   Fires an ActionEvent from the SectionViewer.
   */
  protected void fireEvent(ActionEvent event) {
    // Get the listener list
    listeners2 = actionListeners.getListenerList();
    // Process the listeners last to first
    // List is in pairs, Class and instance
    for (int i = listeners2.length - 2; i >= 0; i -= 2) {
      if (listeners2[i] == ActionListener.class) {
        al = (ActionListener) listeners2[i + 1];
        al.actionPerformed(event);
      }
    }
  } // fireEvent

} // class SectionViewer
