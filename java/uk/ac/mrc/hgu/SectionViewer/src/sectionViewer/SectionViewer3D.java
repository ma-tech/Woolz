package sectionViewer;

import sectionViewer.*;
import hguShape.*;
import hguUntil.*;

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
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.image.codec.jpeg.*;
import com.sun.j3d.utils.image.TextureLoader;

import uk.ac.mrc.hgu.Wlz.*;

public class SectionViewer3D
    extends SectionViewer {

  private final boolean _debug = false;

  threeDMenuHandler handler_7 = null;

  private boolean _3DLoaded = false;
  private float wlzScale = 1.0f;
  private Point3f wlzCentre = null;
  private Point3f wlzFacts = new Point3f(1, 1, 1);
  private Shape3D polygonLine = new Shape3D();
  private Shape3D border = new Shape3D();
  private Arrow arrow = null;
  private boolean viewSolid = false;
  private boolean viewPolygon = false;
  private boolean viewArrow = false;
  private boolean viewTexture = false;
  private String strAxis = "y";
  private Vector3d dir_last = new Vector3d(0.0, -1.0, 0.0);
  private BranchGroup saBG = new BranchGroup();

  private Appearance appearance = new Appearance();
  private LineAttributes lineAttr = new LineAttributes();
  private PolygonAttributes polygenAttr = new PolygonAttributes(
      PolygonAttributes.POLYGON_FILL, PolygonAttributes.CULL_NONE, 0.0f);
  private ColoringAttributes colorAttr = new ColoringAttributes(1.0f, 0.2f,
      0.2f, 0);
  private Color3f[] planeColor = new Color3f[] {gfColor3f.red, gfColor3f.yellow,
      gfColor3f.blue, gfColor3f.pink, gfColor3f.green, gfColor3f.cyan};

  private boolean showFocusSec = false;
  private boolean viewAnatimySV = true;

  //...................
  private String threeDMenuStr = "3D";
  private String threeDMenu_1str = "View Solid";
  private boolean threeDMenu_1state = false;
  private String threeDMenu_2str = "View Arrow";
  private boolean threeDMenu_2state = false;
  private String threeDMenu_3str = "View Border";
  private boolean threeDMenu_3state = true;
  private String threeDMenu_4str = "View Texture";
  private boolean threeDMenu_4state = false;
  private String threeDMenu_5str = "Hide All";
  private String threeDMenu_6str = "Choose Section Color";

  private SVStatus svStatus = new SVStatus();

  private JMenu threeDMenu = new JMenu(threeDMenuStr);
  private JCheckBoxMenuItem threeDMenu_1 = new JCheckBoxMenuItem(
      threeDMenu_1str, threeDMenu_1state);
  private JCheckBoxMenuItem threeDMenu_2 = new JCheckBoxMenuItem(
      threeDMenu_2str, threeDMenu_2state);
  private JCheckBoxMenuItem threeDMenu_3 = new JCheckBoxMenuItem(
      threeDMenu_3str, threeDMenu_3state);
  private JCheckBoxMenuItem threeDMenu_4 = new JCheckBoxMenuItem(
      threeDMenu_4str, threeDMenu_4state);
  private JMenuItem threeDMenu_5 = new JMenuItem(threeDMenu_5str);
  private JMenuItem threeDMenu_6 = new JMenuItem(threeDMenu_6str);

  private PopupMenu pMenu = new PopupMenu();
  private MenuItem menuItemTP = new MenuItem("Add to Tie-Point Table");

  private double meanGreyVal = 1000d;

  //=========================================================
  // constructor
  //=========================================================
  public SectionViewer3D(String viewstr,
                         Container frame,
                         Object parent,
			 BranchGroup obj,
                         boolean showFocusSection) {

    super(viewstr, frame, parent);

    ImageIcon imgIcon = new ImageIcon();
    URL urlIcon =
        SectionViewer3D.class.getResource("images"
                                          + System.getProperty("file.separator")
                                          + "embryo3d.gif");
    if (null != urlIcon){
      imgIcon = new ImageIcon(urlIcon);
      //setIconImage(imgIcon.getImage());
    }

    if (_debug)
      System.out.println("enter SectionViewer3D");

    viewSolid = false;
    viewPolygon = true;

    try {
      Method M1 = _parent.getClass().getDeclaredMethod("threeDLoaded", null);
      _3DLoaded = ( (Boolean) M1.invoke(_parent, null)).booleanValue();
    }
    catch (InvocationTargetException e) {
      System.out.println(e.getMessage());
    }
    catch (NoSuchMethodException ne) {
      System.out.println("threeDLoaded: no such method");
    }
    catch (IllegalAccessException ae) {
      System.out.println(ae.getMessage());
    }

    if (_3DLoaded) {
      try {
        Method M1 = _parent.getClass().getDeclaredMethod("getScale", null);
        wlzScale = ( (Float) M1.invoke(_parent, null)).floatValue();
      }
      catch (InvocationTargetException e) {
        System.out.println(e.getMessage());
      }
      catch (NoSuchMethodException ne) {
        System.out.println("getScale: no such method");
      }
      catch (IllegalAccessException ae) {
        System.out.println(ae.getMessage());
      }

      try {
        Method M1 = _parent.getClass().getDeclaredMethod("getCentre", null);
        wlzCentre = (Point3f) M1.invoke(_parent, null);
      }
      catch (InvocationTargetException e) {
        System.out.println(e.getMessage());
      }
      catch (NoSuchMethodException ne) {
        System.out.println("getCentre: no such method");
      }
      catch (IllegalAccessException ae) {
        System.out.println(ae.getMessage());
      }

      try {
        Method M1 = _parent.getClass().getDeclaredMethod("getFacts", null);
        wlzFacts = (Point3f) M1.invoke(_parent, null);
      }
      catch (InvocationTargetException e) {
        System.out.println(e.getMessage());
      }
      catch (NoSuchMethodException ne) {
        System.out.println("getFacts: no such method");
      }
      catch (IllegalAccessException ae) {
        System.out.println(ae.getMessage());
      }
    }

    saBG.setCapability(BranchGroup.ALLOW_DETACH);

    polygonLine.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
    polygonLine.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
    polygonLine.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);
    polygonLine.setAppearanceOverrideEnable(true);

    border.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
    border.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
    border.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);
    border.setAppearanceOverrideEnable(true);

    lineAttr.setCapability(LineAttributes.ALLOW_WIDTH_WRITE);
    lineAttr.setLineWidth(1.0f);
    lineAttr.setLineAntialiasingEnable(true);

    colorAttr.setCapability(ColoringAttributes.ALLOW_COLOR_WRITE);
    appearance.setCapability(Appearance.ALLOW_COLORING_ATTRIBUTES_WRITE);

    appearance.setLineAttributes(lineAttr);
    appearance.setPolygonAttributes(polygenAttr);
    appearance.setColoringAttributes(colorAttr);
    try {
      colorAttr.setColor(planeColor[(nSV-1) % planeColor.length]);
    }
    catch (Exception exp) {
      colorAttr.setColor(planeColor[0]);
    }
    polygonLine.setAppearance(appearance);

    if(super.colourHandler != null) {
       secColorClt.removeActionListener(super.colourHandler);
    }
    secColorClt.addActionListener(new planeColChooser());

    setShowFocus(showFocusSection);

//...............................
    handler_7 = new threeDMenuHandler();
    threeDMenu.add(threeDMenu_3);
    threeDMenu.add(threeDMenu_1);
    threeDMenu.add(threeDMenu_4);
    threeDMenu.add(threeDMenu_2);
    threeDMenu.add(threeDMenu_5);
    threeDMenu.add(threeDMenu_6);

    threeDMenu_1.addActionListener(handler_7);
    threeDMenu_2.addActionListener(handler_7);
    threeDMenu_3.addActionListener(handler_7);
    threeDMenu_4.addActionListener(handler_7);
    threeDMenu_5.addActionListener(handler_7);
    threeDMenu_6.addActionListener(new planeColChooser());

    _menubar.remove(helpMenu);
    _menubar.add(threeDMenu);
    _menubar.add(helpMenu);

    if (viewType.equals("XY")){
       arrow = new Arrow(new Vector3d(0, 0, 0), new Vector3d(0, 0, 1), "z");
      strAxis = "z";
    } else if (viewType.equals("YZ")){
       arrow = new Arrow(new Vector3d(0, 0, 0), new Vector3d(1, 0, 0), "x");
      strAxis = "x";
    } else if (viewType.equals("ZX")){
      arrow = new Arrow(new Vector3d(0, 0, 0), new Vector3d(0, 1, 0), "y");
      strAxis = "y";
    }
    arrow.setAppearance(appearance);

    saBG.addChild(polygonLine);
    saBG.addChild(border);
    saBG.addChild(arrow);
    obj.addChild(saBG);
    //addWindowListener(new SectionViewer3D_windowAdapter(this));

    if (_debug)
      System.out.println("exit SectionViewer3D");
  } // constructor

//-------------------------------------------------------------
// override methods for opening / closing views
//-------------------------------------------------------------
  public void openView() {
    ViewStructModelTo3DAdaptor WOM23D_1 = null;

    super.openView();

    if (null != WOM23D_1)
      _VSModel.removeChangeListener(WOM23D_1);
    WOM23D_1 = new ViewStructModelTo3DAdaptor(_VSModel, _OBJModel);

    _VSModel.addChangeListener(WOM23D_1);
    /* fireChange is needed to show 3D feedback when view opens */
    _VSModel.fireChange();

    try {
      int[] dstGType = new int[1];
      double[] dstMin = new double[1];
      double[] dstMax = new double[1];
      double[] dstSum = new double[1];
      double[] dstSumSq = new double[1];
      double[] dstMean = new double[1];
      double[] dstStdDev = new double[1];
      int area = WlzObject.WlzGreyStats(_OBJModel.getThreeDObj(),
                                        dstGType, dstMin, dstMax, dstSum,
                                        dstSumSq, dstMean, dstStdDev);
      meanGreyVal = dstMean[0];
    }
    catch (Exception exp) {
      System.out.println(exp.toString());
    }
  } // openView()

//-------------------------------------------------------------
  public void close() {
    super.close();
    saBG.detach();
  }

//---------------------------------------
  // override planeColChooser in SectionViewer
  public class planeColChooser
      implements ActionListener {
    public planeColChooser() {}

    public void actionPerformed(ActionEvent e) {
      Color c = JColorChooser.showDialog(null, "Choose 3D Plane Color",
                                         secColorClt.getBackground());
      if (null != c) {
        secColorClt.setBackground(c);
        colorAttr.setColor(new Color3f(c));
        appearance.setColoringAttributes(colorAttr);
        arrow.setAppearance(appearance);
	updateIntersections();
      }
    }
  }

  public class threeDMenuHandler
      implements ActionListener {
    public threeDMenuHandler() {}

    JCheckBoxMenuItem src;

    public void actionPerformed(ActionEvent e) {
      String valstr;
      if (e.getActionCommand() == threeDMenu_1str) {
        //set view section solid
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != threeDMenu_1state) {
          //set it solid
          viewSolid = true;
          viewPolygon = false;
          viewTexture = false;
          threeDMenu_1.setSelected(true);
          threeDMenu_4.setSelected(false);
          threeDMenu_3.setSelected(false);
        }
        else {
          //set it transparancy
          viewSolid = false;
          threeDMenu_1.setSelected(false);
        }
        try {
          setViewPolygon(_OBJModel, _VSModel);
        }
        catch (Exception exp) {
          exp.printStackTrace();
        }
        return;
      }
      else if (e.getActionCommand() == threeDMenu_2str) {
        //set view direction arrow
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != threeDMenu_2state) {
          //view direction arrow
          viewArrow = true;
          try {
            setViewPolygon(_OBJModel, _VSModel);
          }
          catch (Exception exp) {
            exp.printStackTrace();
          }
          arrow.setWhichChild(Switch.CHILD_ALL);
        }
        else {
          //hide direction arrow
          viewArrow = false;
          arrow.setWhichChild(Switch.CHILD_NONE);
        }
        return;
      }
      else if (e.getActionCommand() == threeDMenu_3str) {
        //set view no plane
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected()) {
          viewSolid = false;
          viewPolygon = true;
          viewTexture = false;
          threeDMenu_1.setSelected(false);
          threeDMenu_4.setSelected(false);
          threeDMenu_3.setSelected(true);

          try {
            setViewPolygon(_OBJModel, _VSModel);
          }
          catch (Exception exp) {
            exp.printStackTrace();
          }
        }
        else {
          viewPolygon = false;
          threeDMenu_3.setSelected(false);
          polygonLine.setGeometry(null);
        }
        return;
      }
      else if (e.getActionCommand() == threeDMenu_4str) {
        //set view no texture
        src = (JCheckBoxMenuItem) e.getSource();
        if (src.isSelected() != threeDMenu_4state) {
          //set texture
          viewSolid = false;
          viewPolygon = false;
          viewTexture = true;
          threeDMenu_1.setSelected(false);
          threeDMenu_4.setSelected(true);
          threeDMenu_3.setSelected(false);
        }
        else {
          //set without texture
          viewTexture = false;
          threeDMenu_4.setSelected(false);
        }
        try {
          setViewPolygon(_OBJModel, _VSModel);
        }
        catch (Exception exp) {
          exp.printStackTrace();
        }
        return;
      }
      else if (e.getActionCommand() == threeDMenu_5str) {
        //Hide all
        threeDMenu_2.setSelected(false);
        threeDMenu_3.setSelected(false);
        threeDMenu_1.setSelected(false);
        threeDMenu_4.setSelected(false);
        viewArrow = false;
        viewSolid = false;
        viewPolygon = false;
        viewTexture = false;
        polygonLine.setGeometry(null);
        border.setGeometry(null);
        arrow.setWhichChild(Switch.CHILD_NONE);
        return;
      }
    } // actionPerformed()
  } // class threeDMenuHandler

//---------------------------------------
// View ADAPTOR
//---------------------------------------
  // change the 3D Polygon when the view structure changes
  public class ViewStructModelTo3DAdaptor
      implements ChangeListener {
    ViewStructModel VSmodel;
    WlzObjModel OBJmodel;

    public ViewStructModelTo3DAdaptor(ViewStructModel mdl1, WlzObjModel mdl2) {
      VSmodel = mdl1;
      OBJmodel = mdl2;
    }

    public void stateChanged(ChangeEvent e) {
      try {
        setViewPolygon(OBJmodel, VSmodel);
      }
      catch (Exception exp) {
        exp.printStackTrace();
      }
    }
  }

//======================================================
  public void set3DLoaded(boolean state) {
    _3DLoaded = state;
  }

  public void setWlzScale(float scale) {
    wlzScale = scale;
  }

  public void setWlzCentre(Point3f Centre) {
    wlzCentre = Centre;
  }

  public void setWlzFacts(Point3f Facts) {
    wlzFacts = Facts;
  }

  void setViewPolygon(WlzObjModel obj, ViewStructModel vs) throws
      Exception {
    WlzObject obj0 = null;
    try {
      obj0 = obj.getThreeDObj();
    }catch( NullPointerException npe){
      return;
    }

    int[] maxVtx = new int[1];
    WlzDVertex3[][] vtxArray = new WlzDVertex3[1][];
    int nVtx = 1;

    //Display polygenlines
    try {
      nVtx = obj0.Wlz3DViewGetBoundingBoxIntersectionA(
          vs.getViewStruct(), maxVtx, vtxArray);
    }
    catch (WlzException exp) {
      exp.printStackTrace();
    }

    if (0 == nVtx)
      return;
    Point3f[] pSection = new Point3f[nVtx];
    int sCount[] = new int[1];
    int i = 0;

    if (showFocusSec)
      lineAttr.setLineWidth(5.0f);
    else
      lineAttr.setLineWidth(3.0f);

    for (i = 0; i < nVtx; i++) {
      pSection[i] = new Point3f(
          (float) (vtxArray[0][i].vtX * wlzFacts.x - wlzCentre.x) / wlzScale,
          (float) (vtxArray[0][i].vtY * wlzFacts.y - wlzCentre.y) / wlzScale,
          (float) (vtxArray[0][i].vtZ * wlzFacts.z - wlzCentre.z) / wlzScale);
    }

    if (viewPolygon) {
      Point3f[] p = new Point3f[pSection.length + 1];
      for (i = 0; i < nVtx; i++)
        p[i] = new Point3f(pSection[i]);

      p[pSection.length] = new Point3f(pSection[0]);

      sCount[0] = p.length;
      LineStripArray line = new LineStripArray(
          p.length, LineStripArray.COORDINATES, sCount);
      line.setCoordinates(0, p);
      polygonLine.setGeometry(line);
      polygonLine.setAppearance(appearance);
      border.setGeometry(null);
    }
    else if (viewSolid) {
      sCount[0] = pSection.length;
      TriangleFanArray triangle = new TriangleFanArray(pSection.length,
          TriangleFanArray.COORDINATES, sCount);
      triangle.setCoordinates(0, pSection);
      polygonLine.setGeometry(triangle);
      polygonLine.setAppearance(appearance);
      border.setGeometry(null);
    }
    else if (viewTexture) {
      try {
        Point3f[] p = new Point3f[pSection.length + 1];
        for (i = 0; i < nVtx; i++)
          p[i] = new Point3f(pSection[i]);
        p[pSection.length] = new Point3f(pSection[0]);
        sCount[0] = p.length;
        LineStripArray line = new LineStripArray(
            p.length, LineStripArray.COORDINATES, sCount);
        line.setCoordinates(0, p);
        border.setAppearance(appearance);
        border.setGeometry(line);

        Appearance texApp = new Appearance();
        PolygonAttributes polyAttr = new PolygonAttributes(
            PolygonAttributes.POLYGON_FILL, PolygonAttributes.CULL_NONE, 0.0f);
        texApp.setPolygonAttributes(polyAttr);

        TextureAttributes textureAttr = new TextureAttributes(
            TextureAttributes.REPLACE, new Transform3D(),
            new Color4f(1.0f, 1.0f, 1.0f, 1.0f), TextureAttributes.NICEST);
        texApp.setTextureAttributes(textureAttr);

        BufferedImage txtImg = null;
        if (viewAnatimySV)
          txtImg = _imgV.getComponentBufferedImage(meanGreyVal);
        else
          txtImg = _imgV.getGreyBufferedImage();

        txtImg = ImageFliperUD(txtImg);
        TextureLoader texLoader = new TextureLoader(txtImg);
        Texture2D texture = (Texture2D) texLoader.getTexture();
        texture.setCapability(Texture.ALLOW_ENABLE_WRITE);
        texture.setBoundaryModeS(Texture.CLAMP);
        texture.setBoundaryModeT(Texture.CLAMP);
        texture.setEnable(true);
        texApp.setTexture(texture);

        WlzThreeDViewStruct VS = _VSModel.getViewStruct();
        double[] xMin = new double[1], yMin = new double[1],
            zMin = new double[1];
        double[] xMax = new double[1], yMax = new double[1],
            zMax = new double[1];
        obj0.Wlz3DViewGetMinvals(VS, xMin, yMin, zMin);
        obj0.Wlz3DViewGetMaxvals(VS, xMax, yMax, zMax);

        float xScale = (float) (xMax[0] - xMin[0]);
        float yScale = (float) (yMax[0] - yMin[0]);

        Point2f texCoord[] = new Point2f[nVtx];
        for (i = 0; i < nVtx; i++) {
          WlzDVertex3[] xyz2D = new WlzDVertex3[1];
          obj0.Wlz3DSectionTransformVtxR(
              vs.getViewStruct(), vtxArray[0][i], xyz2D);

          texCoord[i] = new Point2f( (float) (xyz2D[0].vtX - xMin[0]) / xScale,
                                    (float) (xyz2D[0].vtY - yMin[0]) / yScale);
        }

        sCount[0] = pSection.length;
        TriangleFanArray triangle = new TriangleFanArray(pSection.length,
            TriangleFanArray.COORDINATES |
            TriangleFanArray.TEXTURE_COORDINATE_2,
                          sCount);
        triangle.setCoordinates(0, pSection);
        triangle.setTextureCoordinates(0, texCoord);
        polygonLine.setGeometry(triangle);
        polygonLine.setAppearance(texApp);
      }
      catch (Exception exp) {
        exp.printStackTrace();
      }
    }
    else {
      polygonLine.setGeometry(null);
      border.setGeometry(null);
    }

    //Display arrow
    if (viewArrow) {
      Vector3d p0 = new Vector3d();
      for (i = 0; i < nVtx; i++) {
        p0.x = p0.x + pSection[i].x;
        p0.y = p0.y + pSection[i].y;
        p0.z = p0.z + pSection[i].z;
      }
      p0.x = p0.x / nVtx;
      p0.y = p0.y / nVtx;
      p0.z = p0.z / nVtx;

      Vector3d dir = new Vector3d();
      Vector3d dir1 = new Vector3d(pSection[0].x - pSection[1].x,
                                   pSection[0].y - pSection[1].y,
                                   pSection[0].z - pSection[1].z);
      Vector3d dir2 = new Vector3d(pSection[1].x - pSection[2].x,
                                   pSection[1].y - pSection[2].y,
                                   pSection[1].z - pSection[2].z);
      dir.cross(dir1, dir2);
      if (Math.abs(dir.angle(dir_last)) > 3 * Math.PI / 4)
        dir.negate();

      arrow.setPosition(p0, dir, strAxis);
      dir_last = dir;
    }
  } // end of setView Polygon

  private BufferedImage ImageFliperUD(BufferedImage srcImg) {
    int height = srcImg.getHeight();
    int width = srcImg.getWidth();

    BufferedImage dstImg = new BufferedImage(width, height, srcImg.getType());
    DataBuffer srcBuf = srcImg.getRaster().getDataBuffer();
    DataBuffer dstBuf = dstImg.getRaster().getDataBuffer();

    int srcData;
    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        srcData = srcBuf.getElem(w + (height - h - 1) * width);
        //srcData = (srcData > -600000)?0: srcData;
        dstBuf.setElem(w + h * width, srcData);
      }
    }
    srcImg.flush();
    dstImg.flush();
    return dstImg;
  }

  private BufferedImage ImageTurn180(BufferedImage srcImg) {
    int height = srcImg.getHeight();
    int width = srcImg.getWidth();

    BufferedImage dstImg = new BufferedImage(width, height, srcImg.getType());
    DataBuffer srcBuf = srcImg.getRaster().getDataBuffer();
    DataBuffer dstBuf = dstImg.getRaster().getDataBuffer();

    int srcLen = srcBuf.getSize();
    for (int i = 0; i < srcLen; i++)
      dstBuf.setElem(i, srcBuf.getElem(srcLen - 1 - i));
    srcImg.flush();
    dstImg.flush();

    return dstImg;
  }

  public void setShowFocus(boolean showFocusSec) {
    this.showFocusSec = showFocusSec;
    thisFrameDeactivated();
  }

// override SectionViewer method
  protected void thisFrameActivated() {
    if (showFocusSec) {
      /*
      viewTexture = true;
      viewSolid = false;
      viewPolygon = false;
      threeDMenu_1.setSelected(false);
      */

      try {
        setViewPolygon(_OBJModel, _VSModel);
      }
      catch (Exception exp) {
        exp.printStackTrace();
      }
    }
  }

// override SectionViewer method
  protected void thisFrameDeactivated() {
    if (showFocusSec) {
      /*
      threeDMenu_2.setSelected(false);
      threeDMenu_3.setSelected(false);
      //     threeDMenu_1.setVisible(false);
      //     threeDMenu_4.setVisible(false);
      viewArrow = false;
      */
      polygonLine.setGeometry(null);
      border.setGeometry(null);
    }
  }

  public void showAllSec() {
    viewSolid = svStatus.viewSolid;
    viewPolygon = svStatus.viewPolygon;
    viewTexture = svStatus.viewTexture;
    viewArrow = svStatus.viewArrow;

    if (viewPolygon)
      threeDMenu_3.setSelected(true);
    if (viewSolid)
      threeDMenu_1.setSelected(true);
    if (viewTexture)
      threeDMenu_4.setSelected(true);
    if (viewArrow)
      threeDMenu_2.setSelected(true);

    try {
      setViewPolygon(_OBJModel, _VSModel);
    }
    catch (Exception exp) {
      exp.printStackTrace();
    }

    if (viewArrow) {
      threeDMenu_2.setSelected(true);
      arrow.setWhichChild(Switch.CHILD_ALL);
    }
  }

  public void showFocusSec() {
    svStatus.viewSolid = viewSolid;
    svStatus.viewPolygon = viewPolygon;
    svStatus.viewTexture = viewTexture;
    svStatus.viewArrow = viewArrow;
  }

  public void setAnaChanged(boolean viewAnatimySV) {
    this.viewAnatimySV = viewAnatimySV;
  }

  public void anatomyChange(boolean showInterSecLines) {
    _imgV.setShowInterSecLines(showInterSecLines);
    _imgV.repaint();
    if (viewTexture) {
      try {
        setViewPolygon(_OBJModel, _VSModel);
      }
      catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public void anatomyChange() {
    if (viewTexture) {
      try {
        setViewPolygon(_OBJModel, _VSModel);
      }
      catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public void doIntersection() {
    super.doIntersection();
    this.anatomyChange();
  }

  public void removeIntersection() {
    super.removeIntersection();
    this.anatomyChange();
  }

} // class SectionViewer

/*
class SectionViewer3D_windowAdapter
    extends WindowAdapter {
  SectionViewer3D adaptee;

  SectionViewer3D_windowAdapter(SectionViewer3D adaptee) {
    this.adaptee = adaptee;
  }

  public void windowActivated(WindowEvent e) {
    adaptee.this_windowActivated(e);
  }

  public void windowDeactivated(WindowEvent e) {
    adaptee.this_windowDeactivated(e);
  }
}
*/
