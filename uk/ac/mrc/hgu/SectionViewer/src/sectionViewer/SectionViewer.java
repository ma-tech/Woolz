package sectionViewer;

import sectionViewer.*;
// import hguUntil.*;

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

public class SectionViewer
    extends SectionViewerGUI
    implements Serializable, WlzObjectType, WSetterConstants {

  private final boolean _debug = false;
  private final boolean _NBdebug = false;

  // allows SectionViewer to fire events when true
  private boolean _enabled = true;
  protected boolean _parentIsRoot = true;

  /* parent class must implement the Utils interface */
  protected Container _frame = null;
  protected Object _parent = null;

  protected Cursor xhairCursor = new Cursor(Cursor.CROSSHAIR_CURSOR);
  protected Cursor defCursor = Cursor.getDefaultCursor();

  protected WlzObjModel _OBJModel = null; // grey level
  protected ViewStructModel _VSModel = null;
  public WlzImgView _imgV = null;

  private int _FPBounce = 0; // stops initial false event from PointEntry
  private int _RABounce = 0;

  private WlzObject _anatomyObj = null;
  private WlzObject _constraintObj = null; // wlz obj from thresh constraint

  private AnatomyBuilder _anatBuilder = null;
  private Stack _embFileStack = null;
  private Stack _embObjStack = null;
  private Stack _xembFileStack = null;
  private Stack _xembObjStack = null;
  private Stack _maybeFil = null;
  private Stack _maybeObj = null;

  protected String viewType;
  String _anatomyStr;

  private boolean _anatomy; // for mouse click anatomy
  private boolean _thresholding;
  private boolean _threshConstraint;
  private boolean _fixedPoint = false;
  private boolean _setFixedPoint = false;
  private boolean _showStdCntrls = false;
  private boolean _showUsrCntrls = false;
  private boolean _showIntersection = false;
  private boolean _openingView = false;

  fileMenuHandler handler_1 = null;
  controlMenuHandler handler_2 = null;
  showMenuHandler handler_3 = null;
  helpMenuHandler handler_4 = null;
  // temporarily disabled ... don't remove
  thresholdMenuHandler handler_5 = null;
  //--------------------------------------
  public planeColChooser colourHandler = null;;

  private String SLASH = System.getProperty("file.separator");

  //=========================================================
  // constructor
  //=========================================================
  public SectionViewer(String viewstr, Container frame, Object parent) {

    if (_debug) System.out.println("enter SectionViewer");

    _frame = frame;
    ExtFrameActivationHandler activationHandler1 = new ExtFrameActivationHandler();
    IntFrameActivationHandler activationHandler2 = new IntFrameActivationHandler();

    try {
      Method M1 = _frame.getClass().getDeclaredMethod(
            "setTitle",
            new Class[]{viewstr.getClass( )});
      M1.invoke(_frame, new Object[] {viewstr});
    }
    catch (NoSuchMethodException ne) {
    }
    catch (IllegalAccessException iae) {
    }
    catch (InvocationTargetException ite) {
    }

    try {
      Method M1 = _frame.getClass().getMethod(
            "addWindowListener",
            new Class[]{activationHandler1.getClass().getInterfaces()[0]});
      M1.invoke(_frame, new Object[] {activationHandler1});
    }
    catch (NoSuchMethodException ne) {
    }
    catch (IllegalAccessException iae) {
    }
    catch (InvocationTargetException ite) {
    }

    try {
      Method M1 = _frame.getClass().getDeclaredMethod(
            "addInternalFrameListener",
            new Class[]{activationHandler2.getClass().getInterfaces()[0]});
      M1.invoke(_frame, new Object[] {activationHandler2});
    }
    catch (NoSuchMethodException ne) {
    }
    catch (IllegalAccessException iae) {
    }
    catch (InvocationTargetException ite) {
    }

    _parent = parent;
    try {
      Method M1 = _parent.getClass().getDeclaredMethod("root", null);
    }
    catch (NoSuchMethodException ne) {
      // System.out.println("SVParent sub-classed");
      _parentIsRoot = false;
    }

    _anatomy = false;
    _thresholding = false;
    _threshConstraint = false;

    setCursor(defCursor);
    addFeedback();

/*
    intersectionAdaptor IA = new intersectionAdaptor(this);
    this.addChangeListener(IA);
*/

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
    controlMenu_2_1.addActionListener(handler_2);
    controlMenu_2_2.addActionListener(handler_2);
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

    viewType = viewstr.substring(0, 2);
    if (viewType.equals("XY")) {
      if (pitchSetter != null) {
        pitchSetter.setValue(0.0);
      }
      if (yawSetter != null) {
        yawSetter.setValue(0.0);
      }
      if (rollSetter != null) {
        if (rollSetter.isSliderEnabled() == true) {
          rollSetter.setValue(0.0);
        }
        else {
          rollSetter.setSliderEnabled(true);
          rollSetter.setValue(0.0);
          rollSetter.setSliderEnabled(false);
        }
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
        if (rollSetter.isSliderEnabled() == true) {
          rollSetter.setValue(90.0);
        }
        else {
          rollSetter.setSliderEnabled(true);
          rollSetter.setValue(90.0);
          rollSetter.setSliderEnabled(false);
        }
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
        if (rollSetter.isSliderEnabled() == true) {
          rollSetter.setValue(90.0);
        }
        else {
          rollSetter.setSliderEnabled(true);
          rollSetter.setValue(90.0);
          rollSetter.setSliderEnabled(false);
        }
      }
    }
    resetFeedbackText();


    if (_debug)
      System.out.println("exit SectionViewer");
  } // constructor

//-------------------------------------------------------------
// methods for adding / removing panels from the gui
//-------------------------------------------------------------
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
    //_rootPane.validate();
    revalidate();
    repaint();
  }

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
    //_rootPane.validate();
    revalidate();
    repaint();
  }

  //...............................

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
    //_rootPane.validate();
    revalidate();
    repaint();
  }

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
    //_rootPane.validate();
    revalidate();
    repaint();
  }

  //...............................

  protected void addFeedback() {
    feedbackImagePanel.add(feedbackPanel, BorderLayout.NORTH);
    //_rootPane.validate();
    revalidate();
    repaint();
  }

  protected void removeFeedback() {
    feedbackImagePanel.remove(feedbackPanel);
    feedbackImagePanel.repaint();
    //_rootPane.validate();
    revalidate();
    repaint();
  }

//-------------------------------------------------------------
// methods for opening / closing views
//-------------------------------------------------------------
  public void openView() {
    if (null == _parent) {
      System.out.println("There no object in system!!!");
      return;
    }

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
      if (_parentIsRoot) {
        M1 = _parent.getClass().getDeclaredMethod("getOBJModel", null);
      }
      else {
        M1 = _parent.getClass().getSuperclass().getDeclaredMethod("getOBJModel", null);
      }
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
    try {
      Method M1 = null;
      if (_parentIsRoot) {
        M1 = _parent.getClass().getDeclaredMethod("getTitleText", null);
      }
      else {
        M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
            "getTitleText", null);
      }
      setTitleText( (String) M1.invoke(_parent, null));
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("getTitleText: no such method 1");
      System.out.println(ne.getMessage());
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }
    //setTitleText(_parent.getTitleText());
    setViewTitle();
    try {
      Method M1 = null;
      if (_parentIsRoot) {
        M1 = _parent.getClass().getDeclaredMethod("getAnatomyBuilder", null);
      }
      else {
        M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
            "getAnatomyBuilder", null);
      }
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
  protected void setDistLimits(double val) {

    double results[];
    double min = 0.0;
    double max = 0.0;

    results = _OBJModel.getMaxMin(_VSModel.getViewStruct());
    min = results[5];
    max = results[2];

    distSetter.setMin(min);
    distSetter.setMax(max);
    distSetter.setValue(val);
  }

//-------------------------------------------------------------
// set the title text
//-------------------------------------------------------------
  protected void setTitleText(String title) {
    titleText.setText(title);
  }

//-------------------------------------------------------------
// re-set Feedback text
//-------------------------------------------------------------
  protected void resetFeedbackText() {
    resetPosGreyText();
    resetAnatomyText();
  }

//-------------------------------------------------------------
// re-set pos & grey text
//-------------------------------------------------------------
  protected void resetPosGreyText() {
    xyzTextField.setText("-----");
    valueTextField.setText("-");
  }

//-------------------------------------------------------------
// re-set anatomy text
//-------------------------------------------------------------
  protected void resetAnatomyText() {
    anatomyTextField.setText("-----");
  }

//-------------------------------------------------------------
// print model info
//-------------------------------------------------------------
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
  protected boolean isValidName(String dir, String fil) {
    if (fil.endsWith(".wlz")) {
      return true;
    }
    else {
      return false;
    }
  }

  //-------------------------------------------------------------
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
  public void anatomyFromMenu() {

    AnatomyElement anatArr[] = null;
    if (null != _parent) {
      try {
        Method M1 = null;
        if (_parentIsRoot) {
          M1 = _parent.getClass().getDeclaredMethod("getAnatomyArr", null);
        }
        else {
          M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
              "getAnatomyArr", null);
        }
        anatArr = (AnatomyElement[]) M1.invoke(_parent, null);
      }
      catch (InvocationTargetException e) {
        System.out.println(e.getMessage());
      }
      catch (NoSuchMethodException ne) {
        System.out.println("getAnatomyArray: no such method 1");
        System.out.println(ne.getMessage());
      }
      catch (IllegalAccessException ae) {
        System.out.println(ae.getMessage());
      }
      anatomyFromMenu(anatArr);
    }
  }

  public void anatomyFromMenu(AnatomyElement[] arr) {
    int num = AnatKey._nrows;
    WlzObject obj2D[] = new WlzObject[num];
    boolean viz[] = new boolean[num];

    for (int i = 0; i < num; i++) {
      if ( (arr[i] != null) && (arr[i].isRemoved() == false)) {
        viz[i] = arr[i].isVisible();
        obj2D[i] = _OBJModel.makeSection(
            arr[i].getObj(),
            _VSModel.getViewStruct());
        if (obj2D[i] != null) {
          if (WlzObject.WlzObjGetType(obj2D[i]) ==
              WlzObjectType.WLZ_EMPTY_OBJ) {
            obj2D[i] = null;
          }
        }
      }
      else {
        obj2D[i] = null;
        viz[i] = false;
      }
    }

    try {
      _imgV.clearAnatomy();
      _imgV.setAnatomyObj(obj2D, viz);
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
  protected Line2D.Double[] getIntersectionArr() {

    Line2D.Double intersectionArr[] = null;
    SectionViewer SV = null;
    Vector svVec = null;
    try {
      Method M1 = null;
      if (_parentIsRoot) {
        M1 = _parent.getClass().getDeclaredMethod("getOpenViews", null);
      }
      else {
        M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
            "getOpenViews", null);
      }
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

      //printMaxMin(xmin[0],ymin[0],xmax[0],ymax[0]);

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
          xp = vtx.vtX;
          yp = vtx.vtY;
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
	  j++;
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
  protected Color[] getInterColArr() {

    Color colarr[] = null;
    SectionViewer SV = null;
    Vector svVec = null;
    try {
      Method M1 = null;
      if (_parentIsRoot) {
        M1 = _parent.getClass().getDeclaredMethod("getOpenViews", null);
      }
      else {
        M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
            "getOpenViews", null);
      }
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
  protected void doIntersection() {
    if (_imgV == null) return;
    if (_showIntersection) {
      _imgV.setIntersectionVec(getIntersectionArr());
      _imgV.setInterColVec(getInterColArr());
    }
    _imgV.repaint();
  }

//-------------------------
  protected void updateIntersections() {

    SectionViewer SV = null;
    Vector svVec = null;
    try {
      Method M1 = null;
      if (_parentIsRoot) {
        M1 = _parent.getClass().getDeclaredMethod("getOpenViews", null);
      }
      else {
        M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
            "getOpenViews", null);
      }
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
  protected void printIntersection(Line2D.Double line) {

     String x1Str;
     String y1Str;
     String x2Str;
     String y2Str;

     if(line == null) return;

     x1Str = Double.toString(line.getX1());
     x2Str = Double.toString(line.getX2());
     y1Str = Double.toString(line.getY1());
     y2Str = Double.toString(line.getY2());

     System.out.print("intersection line: "+x1Str+","+y1Str);
     System.out.println(" --- "+x2Str+","+y2Str);

  }

//-------------------------
  protected void printMaxMin(double x1, double y1, double x2, double y2) {

     String x1Str;
     String y1Str;
     String x2Str;
     String y2Str;

     x1Str = Double.toString(x1);
     y1Str = Double.toString(y1);
     x2Str = Double.toString(x2);
     y2Str = Double.toString(y2);

     System.out.print("x range: "+x1+","+x2);
     System.out.println("y range: "+y1+","+y2);

  }

//-------------------------
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

//-------------------------------------------------------------
// remove overlay
//-------------------------------------------------------------
  protected void removeOverlay() {
    _imgV.clearOverlay();
    _maybeObj = null;
    _maybeFil = null;
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove intersection
//-------------------------------------------------------------
  protected void removeIntersection() {
    _imgV.clearIntersection();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove threshold
//-------------------------------------------------------------
  protected void removeThreshold() {
    _imgV.clearThreshold();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove anatomy
//-------------------------------------------------------------
  public void removeAnatomy() {
    _imgV.clearAnatomy();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove threshConstraint
//-------------------------------------------------------------
  protected void removeThreshConstraint() {
    _imgV.clearThreshConstraint();
    _imgV.repaint();
  }

//-------------------------------------------------------------
// remove fixedPoint
//-------------------------------------------------------------
  protected void removeFixedPoint() {
    _imgV.enableFixedPoint(false);
    _imgV.repaint();
  }

//-------------------------------------------------------------
// set view title
//-------------------------------------------------------------
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
    String viewTitle = pitch + " | " + yaw + " | " + roll;
    try {
      Method M1 = _frame.getClass().getMethod(
            "setTitle",
            new Class[]{viewTitle.getClass( )});
      M1.invoke(_frame, new Object[] {viewTitle});
    }
    catch (NoSuchMethodException ne) {
    }
    catch (IllegalAccessException iae) {
    }
    catch (InvocationTargetException ite) {
    }

  }

//-------------------------------------------------------------
// save image as jpeg
//-------------------------------------------------------------
  protected void saveImage(String str) {
    File fil = new File(str);
    BufferedImage img = null;

    img = _imgV.getComponentBufferedImage(1000d);

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
  protected double zf2mag(double zf) {
    //return (zf < 0.0) ? 1.0 + 0.001 * zf : 1.0 + 0.01 * zf;
    return zf/100.0;
  }

//-------------------------------------------------------------
// checkFilename for extension
//-------------------------------------------------------------
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
  public void setViewMode(String mode) {
    if (_VSModel != null) {
      _VSModel.setViewMode(mode);
      if (mode.equals("absolute") == true) {
        controlMenu_2_2.setSelected(true);
        rollSetter.setSliderEnabled(true);
      }
      else {
        controlMenu_2_1.setSelected(true);
        rollSetter.setSliderEnabled(false);
      }
    }
  }

//-------------------------------------------------------------
// get view structure model etc
//-------------------------------------------------------------
  protected ViewStructModel getViewStructModel() {
    return _VSModel;
  }

  protected WlzObjModel getWlzObjModel() {
    return _OBJModel;
  }

  protected WlzImgView getImageView() {
    return _imgV;
  }

  protected Zoom getZoomSetter() {
    return zoomSetter;
  }

  protected WSetter getDistSetter() {
    return distSetter;
  }

  protected WSetter getYawSetter() {
    return yawSetter;
  }

  protected WSetter getPitchSetter() {
    return pitchSetter;
  }

  protected WSetter getRollSetter() {
    return rollSetter;
  }

  protected JButton getSecColorClt() {
    return secColorClt;
  }

//-------------------------------------------------------------
// enable / disable thresholding & threshConstraint
//-------------------------------------------------------------
  protected void enableThresholding(boolean state) {
    _thresholding = state;
  }

//-------------------------------------------------------------
  protected void enableThreshConstraint(boolean state) {
    _threshConstraint = state;
    if (state) {
      if (!thresholdMenu_3state) {
        thresholdMenu_3.doClick();
      }
    }
  }

//-------------------------------------------------------------
  protected void disableAnatomy() {
    _anatomy = false;
  }

  protected void enableAnatomy() {
    _anatomy = true;
  }

//-------------------------------------------------------------
// set anatomy text field for menu selected anatomy
//-------------------------------------------------------------
  protected void setAnatomyText(String str) {
    resetAnatomyText();
    anatomyTextField.setText(str);
  }

//-------------------------------------------------------------
// reset fixed point
//-------------------------------------------------------------
  protected void resetFixedPoint(boolean showFP) {

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
// reset user interface controls
//-------------------------------------------------------------
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
      if (rollSetter.isSliderEnabled()) {
        rollSetter.setValue(0.0);
      }
      else {
        rollSetter.setSliderEnabled(true);
        rollSetter.setValue(0.0);
        rollSetter.setSliderEnabled(false);
      }
    }
    else if (viewType.equals("YZ")) {
      pitchSetter.setValue(90.0);
      yawSetter.setValue(0.0);
      if (rollSetter.isSliderEnabled()) {
        rollSetter.setValue(90.0);
      }
      else {
        rollSetter.setSliderEnabled(true);
        rollSetter.setValue(90.0);
        rollSetter.setSliderEnabled(false);
      }
    }
    else if (viewType.equals("ZX")) {
      pitchSetter.setValue(90.0);
      yawSetter.setValue(90.0);
      if (rollSetter.isSliderEnabled()) {
        rollSetter.setValue(90.0);
      }
      else {
        rollSetter.setSliderEnabled(true);
        rollSetter.setValue(90.0);
        rollSetter.setSliderEnabled(false);
      }
    }
    else {
      System.out.println("unknown view type in resetGUI(): " + viewType);
    }
    distSetter.setValue(0.0);
    zoomSetter.setValue(100);
    setViewMode("up_is_up");
    removeThreshold();
    removeThreshConstraint();
    enableThresholding(false);
    enableThreshConstraint(false);
    disableAnatomy();
    resetFeedbackText();
    resetFixedPoint(false);
    setSVEnabled(true);

    if (_debug)
      System.out.println("exit resetGUI");
  }

//-------------------------------------------------------------
// constraint to Wlz
//-------------------------------------------------------------
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
  protected void close() {
    Object params[] = {this};
    Vector openViews = null;
    try {
      Method M1 = null;
      if (_parentIsRoot) {
        M1 = _parent.getClass().getDeclaredMethod(
	    "removeView",
            new Class[]{this.getClass( )});
      }
      else {
        M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
	    "removeView",
            new Class[]{this.getClass( ).getSuperclass()});
      }
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
      Method M1 = _frame.getClass().getMethod(
            "dispose", null);
      M1.invoke(_frame, null);
    }
    catch (NoSuchMethodException ne) {
    }
    catch (IllegalAccessException iae) {
    }
    catch (InvocationTargetException ite) {
    }
  }

//-------------------------------------------------------------
// deal with container activation de-activation
//-------------------------------------------------------------
  protected void thisFrameActivated() {
  }

  protected void thisFrameDeactivated() {
  }

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//-------------------------------------------------------------
// connect up the adaptors
//-------------------------------------------------------------
  ZoomToZoomAdaptor Z2Z_1;
  WSetterToDistAdaptor W2D_1;
  WSetterToPitchAdaptor W2P_1;
  WSetterToYawAdaptor W2Y_1;
  WSetterToRollAdaptor W2R_1;
  ViewStructModelToViewAdaptor WOM2V_1;
  WlzImgViewToPosAdaptor WI2P_1;
  WlzImgViewToGreyAdaptor WI2G_1;
  WlzImgViewToAnatAdaptor WI2A_1;
  WlzImgViewToFPAdaptor WI2FP_1;
  MouseToFBAdaptor M2FB_1;
  MouseToThresholdAdaptor M2T_1;
  MouseToThreshConstraintAdaptor M2TC_1;

  // adaptors for grey level
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
    if (WI2FP_1 != null)
      _imgV.removeChangeListener(WI2FP_1);
    WI2FP_1 = new WlzImgViewToFPAdaptor(_imgV, _OBJModel, _VSModel);
    _imgV.addChangeListener(WI2FP_1);
    //...............................
    if (M2FB_1 != null) {
      _imgV.removeMouseListener(M2FB_1);
      _imgV.removeMouseMotionListener(M2FB_1);
    }
    M2FB_1 = new MouseToFBAdaptor(_imgV);
    _imgV.addMouseListener(M2FB_1);
    _imgV.addMouseMotionListener(M2FB_1);
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
    //--------------------------------------

    setSVEnabled(true);
    if (_debug)
      System.out.println("exit connectAdaptors_1");
  } //connectAdaptors_1()

  protected void setSVEnabled(boolean state) {
    _enabled = state;
  }

  protected boolean isSVEnabled() {
    return _enabled;
  }

//==========================================  
// get and set methods
//==========================================  

   public double getDist() {
      double ret = 0;
      if(distSetter != null) {
         ret = ((Double)distSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   public void setDist(double val) {
      if(distSetter != null) {
         distSetter.setValue(val);
      }
   }
//......................................
   public double getPitch() {
      double ret = 0;
      if(pitchSetter != null) {
         ret = ((Double)pitchSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   public void setPitch(double val) {
      if(pitchSetter != null) {
         pitchSetter.setValue(val);
      }
   }
//......................................
   public double getYaw() {
      double ret = 0;
      if(yawSetter != null) {
         ret = ((Double)yawSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   public void setYaw(double val) {
      if(yawSetter != null) {
         yawSetter.setValue(val);
      }
   }
//......................................
   public double getRoll() {
      double ret = 0;
      if(rollSetter != null) {
         ret = ((Double)rollSetter.getValue().elementAt(0)).doubleValue();
      }
      return ret;
   }

   public void setRoll(double val) {
      /* roll should not be adjusted if view mode is 'up_is_up' */
      if((rollSetter != null)&&(controlMenu_2_2.isSelected())) {
         rollSetter.setValue(val);
      }
   }
//......................................
   public int getZoom() {
      int ret = 0;
      if(zoomSetter != null) {
         ret = zoomSetter.getValue();
      }
      return ret;
   }

   public void setZoom(int val) {
      if(zoomSetter != null) {
         zoomSetter.setValue(val);
      }
   }
//......................................

//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
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
          if (src.isSelected() == true) {
            addStandardControls();
            _showStdCntrls = true;
          }
          else {
            removeStandardControls();
            _showStdCntrls = false;
          }
          controlMenu_1_1state = src.isSelected();
        }
      }
      else if (e.getActionCommand().equals(controlMenu_1_2str)) {
        // see user defined rotation controls
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != controlMenu_1_2state) {
          if (src.isSelected() == true) {
            addUserControls();
            _showUsrCntrls = true;
          }
          else {
            removeUserControls();
            _showUsrCntrls = false;
          }
          controlMenu_1_2state = src.isSelected();
        }
      }
      // view mode ------------------------------
      else if (e.getActionCommand().equals(controlMenu_2_1str)) {
        // up is up mode => roll is fixed
        setViewMode(controlMenu_2_1str);
        rollSetter.setSliderEnabled(true);
        rollSetter.setValue(0.0);
        rollSetter.setSliderEnabled(false);
      }
      else if (e.getActionCommand().equals(controlMenu_2_2str)) {
        setViewMode(controlMenu_2_2str);
      }
      // fixed point ------------------------------
      else if (e.getActionCommand().equals(controlMenu_3_1str)) {
        // change using mouse
        _setFixedPoint = true;
        setCursor(xhairCursor);
      }
      else if (e.getActionCommand().equals(controlMenu_3_2str)) {
        // change by typing coordinates
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
      // rotation axis ------------------------------
      else if (e.getActionCommand().equals(controlMenu_4_2str)) {
        // change by typing coordinates
        _VSModel.getFixedPoint(xa, ya, za);
        vals = new double[3];
        vals[0] = xa[0];
        vals[1] = ya[0];
        vals[2] = za[0];
        _RABounce = 0;
        pe = new PointEntry("Rotation Axis Coordinates");
        pe.setSize(250, 250);
        pe.pack();
        pe.setVisible(true);
        pe.setInitVals(vals);
        pe.addChangeListener(new RotationAxisAdaptor(pe));
        vals = null;
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
      // rotation axis ------------------------------
      else if (e.getActionCommand().equals(showMenu_5str)) {
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != showMenu_5state) {
          if (src.isSelected()) {
          }
          else {
          }
          showMenu_5state = src.isSelected();
        }
      }
    }
  } // class showMenuHandler
//---------------------------------------
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
          if (_parentIsRoot) {
            M1 = _parent.getClass().getDeclaredMethod("getHelpBroker", null);
          }
          else {
            M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
                "getHelpBroker", null);
          }
          hb = (HelpBroker) M1.invoke(_parent, null);
        }
        catch (InvocationTargetException ie) {
          System.out.println(ie.getMessage());
        }
        catch (NoSuchMethodException ne) {
          System.out.println("getHelpBroker: no such method");
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
  public class thresholdMenuHandler
      implements ActionListener {

    /* the options are not part of a button group as they
       are not mutually exclusive */
    JCheckBoxMenuItem src;
    JOptionPane jop;
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

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//---------------------------------------
// Control ADAPTORS, wsetter to wlz object
//---------------------------------------
  // change 'distance' using WSetter control
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
    }
  }

//---------------------------------------
  // change 'pitch' using WSetter control
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
      if (!control.isAdjusting()) {
        VSmodel.getDist(valArr);
        setDistLimits(valArr[0]);
      }
    }
  }

//---------------------------------------
  // change 'yaw' using WSetter control
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
// View ADAPTORS
//---------------------------------------
  // change the displayed section when the view structure changes
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
    AnatomyElement arr[];
    WlzObject obj2D[];

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

      int num = AnatKey._nrows;

      obj2D = new WlzObject[num];
      boolean viz[] = new boolean[num];

      for (int i = 0; i < num; i++) {
        obj2D[i] = null;
        viz[i] = false;
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
//-------------------------
      // draw the anatomy components
      try {
        if (null != _parent) {
          try {
            Method M1 = null;
            if (_parentIsRoot) {
              M1 = _parent.getClass().getDeclaredMethod("getAnatomyArr", null);
            }
            else {
              M1 = _parent.getClass().getSuperclass().getDeclaredMethod(
                  "getAnatomyArr", null);
            }
            arr = (AnatomyElement[]) M1.invoke(_parent, null);
          }
          catch (InvocationTargetException ie) {
            System.out.println(ie.getMessage());
          }
          catch (NoSuchMethodException ne) {
            System.out.println("getAnatomyArr: no such method");
            System.out.println(ne.getMessage());
          }
          catch (IllegalAccessException ae) {
            System.out.println(ae.getMessage());
          }
        }

        if ( (arr != null) && (viz != null) && (obj2D != null)) {
          for (int i = 0; i < num; i++) {
            if (arr[i] != null) {
              viz[i] = (arr[i].isVisible() &&
                        !arr[i].isRemoved());
              obj2D[i] = _OBJModel.makeSection(
                  arr[i].getObj(),
                  VS);
              if (obj2D[i] != null) {
                if (WlzObject.WlzObjGetType(obj2D[i]) ==
                    WlzObjectType.WLZ_EMPTY_OBJ) {
                  obj2D[i] = null;
                }
              }
            }
            else {
              viz[i] = false;
              obj2D[i] = null;
            }
          }
          view.setAnatomyObj(obj2D, viz);
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
      //_rootPane.validate();
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
  // connect mouse events to setting of fixed point
  public class WlzImgViewToFPAdaptor
      implements ChangeListener {

    WlzImgView ImgModel;
    WlzObjModel ObjModel;
    ViewStructModel VSModel;
    WlzThreeDViewStruct VS;
    Point pos;
    double pt3d[];

    public WlzImgViewToFPAdaptor(WlzImgView mdl1,
                                 WlzObjModel mdl2,
                                 ViewStructModel mdl3) {
      ImgModel = mdl1;
      ObjModel = mdl2;
      VSModel = mdl3;
      VS = VSModel.getViewStruct();
    }

    public void stateChanged(ChangeEvent e) {
      if (!_setFixedPoint)
        return;
      pos = ImgModel.getPos();
      pt3d = ObjModel.get3DPoint(pos, VS);
      VSModel.setFixedPoint(pt3d[0],
                            pt3d[1],
                            pt3d[2]);
      setDistLimits(0.0);
      _setFixedPoint = false;
      setCursor(defCursor);
      Runnable fpClick = new Runnable() {
        public void run() {
          if (!showMenu_4state)
            showMenu_4.doClick();
        }
      };
      SwingUtilities.invokeLater(fpClick);
    }

  }

//---------------------------------------
  // change 'zoom' using Zoom control
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
          view.setThreshConstraintOffsets(bbox2.xMin,
                                          bbox2.yMin);
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
  // connects window events to focus notification
  public class WindowFocusHandler
      implements WindowFocusListener {

    ActionEvent event = null;
    String actionCommand = "focus";

    public WindowFocusHandler() {
    }

    public void windowLostFocus(WindowEvent e) {
    }

    public void windowGainedFocus(WindowEvent e) {
      event = new ActionEvent(e.getSource(),
                              ActionEvent.ACTION_PERFORMED,
                              actionCommand);
      fireEvent(event);
    }
  }

//---------------------------------------
  // notifies window activation events
  public class ExtFrameActivationHandler
      implements WindowListener {

    public ExtFrameActivationHandler() {
    }

    public void windowActivated(WindowEvent e) {
       thisFrameActivated();
    }

    public void windowDeactivated(WindowEvent e) {
       thisFrameDeactivated();
    }

    public void windowOpened(WindowEvent e) { }
    public void windowClosing(WindowEvent e) { }
    public void windowClosed(WindowEvent e) { }
    public void windowIconified(WindowEvent e) { }
    public void windowDeiconified(WindowEvent e) { }

  }

//---------------------------------------
  // notifies internal frame activation events
  public class IntFrameActivationHandler
      implements InternalFrameListener {

    public IntFrameActivationHandler() {
    }

    public void internalFrameActivated(InternalFrameEvent e) {
       thisFrameActivated();
    }

    public void internalFrameDeactivated(InternalFrameEvent e) {
       thisFrameDeactivated();
    }

    public void internalFrameOpened(InternalFrameEvent e) { }
    public void internalFrameClosing(InternalFrameEvent e) { }
    public void internalFrameClosed(InternalFrameEvent e) { }
    public void internalFrameIconified(InternalFrameEvent e) { }
    public void internalFrameDeiconified(InternalFrameEvent e) { }

  }

//-------------------------------------------------------------
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
  public class RotationAxisAdaptor
      implements ChangeListener {

    double vals[];
    PointEntry _pe;

    public RotationAxisAdaptor(PointEntry pe) {
      _pe = pe;
    }

    public void stateChanged(ChangeEvent e) {
    }
  }

//-------------------------------------------------------------

//=============================================================
//-------------------------------------------------------------
// thread stuff
//-------------------------------------------------------------

  Object _lock = new Object();
//-------------------------------------------------------------
  // stack
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
  protected EventListenerList changeListeners =
      new EventListenerList();

//-------------------------------------------------------------
  // add a listener to the register
  public void addChangeListener(ChangeListener x) {
    changeListeners.add(ChangeListener.class, x);

    // bring it up to date with current state
    x.stateChanged(new ChangeEvent(this));
  }

//-------------------------------------------------------------
  // remove a listener from the register
  public void removeChangeListener(ChangeListener x) {
    changeListeners.remove(ChangeListener.class, x);
  }

//-------------------------------------------------------------
  private ChangeEvent ce;
  private Object[] listeners;
  private ChangeListener cl;
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
  protected EventListenerList actionListeners =
      new EventListenerList();

//-------------------------------------------------------------
  // add a listener to the register
  public void addActionListener(ActionListener x) {
    actionListeners.add(ActionListener.class, x);
  }

//-------------------------------------------------------------
  // remove a listener from the register
  public void removeActionListener(ActionListener x) {
    actionListeners.remove(ActionListener.class, x);
  }

//-------------------------------------------------------------
  private Object[] listeners2;
  private ActionListener al;
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
