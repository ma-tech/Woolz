package sectionViewer;

import sectionViewer.*;

import java.io.*;
import java.util.*;
import javax.vecmath.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.media.j3d.*;

public class SVParent3D
    extends SVParent {

  private final boolean debug = false;

  private boolean _3DLoaded = false;
  private float _scale = 1.0f;
  private Point3f _centre = null;
  private Point3f _facts = null;

  /**
   * Purpose:    Constructor
   **/
  public SVParent3D() {
    super();
    init2D();
    initHelp();
  }

  //-------------------------------------------------------------
  protected void init2D() {
    super.init2D();

    if (K2C_1 != null)
      _key.removeActionListener(K2C_1);
    K2C_1 = new keyToControllerAdaptor3D();
    _key.addActionListener(K2C_1);
  }

  //-------------------------------------------------------------
  public JFrame createExternalView(String viewStr,
				  BranchGroup obj,
				  boolean showFocusSection) {

     return createExternalView(viewStr,
                               obj,
		               showFocusSection,
			       true);
  }
  //-------------------------------------------------------------
  public JFrame createExternalView(String viewStr,
				  BranchGroup obj,
				  boolean showFocusSection,
				  boolean visible) {

    JFrame extFr = null;
    JPanel pan1 = null;
    SectionViewer3D SV = null;

    pan1 = new JPanel();
    pan1.setLayout(new BorderLayout());

    extFr = new JFrame("view");
    extFr.setContentPane(pan1);

    SV = new SectionViewer3D(viewStr,
                             extFr,
			     this,
			     obj,
			     showFocusSection);
    if(SV != null) {
       addView(SV);
       SV.openView();
       pan1.add(SV, BorderLayout.CENTER);
       extFr.pack();
       extFr.setVisible(visible);
    }

    return extFr;
  }

  //-------------------------------------------------------------
  public JInternalFrame createInternalView(String viewStr,
				           BranchGroup obj,
				           boolean showFocusSection) {

     return createInternalView(viewStr,
                               obj,
			       defViewW,
			       defViewH,
		               showFocusSection,
			       true);
  }
  //-------------------------------------------------------------
  public JInternalFrame createInternalView(String viewStr,
			                   BranchGroup obj,
					   int W,
					   int H,
				           boolean showFocusSection,
				           boolean visible) {

    JInternalFrame intFr = null;
    JPanel pan1 = null;
    SectionViewer3D SV = null;

    pan1 = new JPanel();
    pan1.setLayout(new BorderLayout());

    intFr = new JInternalFrame("view", true, true, true, true);
    intFr.setContentPane(pan1);

    SV = new SectionViewer3D(viewStr,
                             intFr,
			     this,
			     obj,
			     showFocusSection);
    if(SV != null) {
       addView(SV);
       SV.openView();
       pan1.add(SV, BorderLayout.CENTER);
       //intFr.pack();
       intFr.setPreferredSize(new Dimension(W, H));
       intFr.setVisible(visible);
    }

    return intFr;
  }

  //-------------------------------------------------------------
  public void set3DParams(float scale, Point3f centre, Point3f facts) {
    if ( (centre == null) || (facts == null)) {
      System.out.println("returning early");
      return;
    }
    _scale = scale;
    _centre = centre;
    _facts = facts;
    _3DLoaded = true;

    /* if there are any open views, sort them out */
    SectionViewer3D SV = null;
    int numViews = _openViews.size();
    if (numViews <= 0)
      return;

    for (int i = 0; i < numViews; i++) {
      //System.out.println("set3DParams for view "+i);
      SV = (SectionViewer3D) _openViews.elementAt(i);
      SV.setWlzScale(scale);
      SV.setWlzCentre(centre);
      SV.setWlzFacts(facts);
      SV.set3DLoaded(true);
    }
  }

  //-------------------------------------------------------------
  public boolean threeDLoaded() {
    /* for views opened after 'set3DParams was called */
    return _3DLoaded;
  }

  //-------------------------------------------------------------
  public float getScale() {
    /* for views opened after 'set3DParams was called */
    return _scale;
  }

  //-------------------------------------------------------------
  public Point3f getCentre() {
    /* for views opened after 'set3DParams was called */
    return _centre;
  }

  //-------------------------------------------------------------
  public Point3f getFacts() {
    /* for views opened after 'set3DParams was called */
    return _facts;
  }

  //-------------------------------------------------------------
  /*
     if this class is sub-classed again it will be necessary to provide
     the following methods so that reflection from SectionViewer can
     find them...
     getOBJModel()
     getTitleText()
     getAnatomyBuilder()
     getAnatomyArr()
     getOpenViews()
     getHelpBroker()
   */
  public class keyToControllerAdaptor3D
      extends keyToControllerAdaptor
      implements ActionListener {

    int numViews;
    int indx;
    SectionViewer3D SV = null;
    Color newCol = null;

    public keyToControllerAdaptor3D() {
    }

    public void actionPerformed(ActionEvent e) {
      String cmd = e.getActionCommand();
      boolean wasRemoved = false;

      indx = (
          Integer.valueOf(cmd.substring(cmd.length() - 1))).intValue();
      //if(_anatomyArr[indx] == null) return;

      if (cmd.indexOf("makeVisible") != -1) {
        if (_anatomyArr[indx] != null) {
          _anatomyArr[indx].setVisible(true);
        }
      }
      else if (cmd.indexOf("makeInvisible") != -1) {
        if (_anatomyArr[indx] != null)
          _anatomyArr[indx].setVisible(false);
      }
      else if (cmd.indexOf("zap") != -1) {
        if (_anatomyArr[indx] != null) {
          wasRemoved = _anatomyArr[indx].isRemoved();
          _key.setZapToolTips(indx, !wasRemoved);
          _anatomyArr[indx].setRemoved(!wasRemoved);
        }
      }
      else if (cmd.indexOf("colour") != -1) {
        newCol = JColorChooser.showDialog(null,
                                          "choose colour for anatomy component",
                                          new Color(_key._cols[indx]));
        _key.setCol(newCol, indx);
      }
      else {
        System.out.println("unknown command");
      }

      AnatKey.update(_anatomyArr);
      doChange();
    }

    private void doChange() {
      numViews = _openViews.size();
      if (numViews == 0)
        return;

      for (int i = 0; i < numViews; i++) {
        SV = (SectionViewer3D) _openViews.elementAt(i);
        SV.removeAnatomy();
        SV.anatomyFromMenu(_anatomyArr);
        SV.anatomyChange();
      }
    }
  }
}
