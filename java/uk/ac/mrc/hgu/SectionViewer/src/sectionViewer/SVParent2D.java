package sectionViewer;
import sectionViewer.*;

import java.net.*;
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.JColorChooser;
import javax.help.*;

import uk.ac.mrc.hgu.Wlz.*;

public class SVParent2D implements SVParent {

  private String _titleText = "";
  private File selectedGreyFile = null;

  private WlzObjModel _OBJModel = null;
  private AnatomyBuilder _anatBuilder = null;
  private SVLocks _Locks = null;

  protected Vector _openViews = null;

  private Stack _embTreeStack = null;
  private Stack _xembTreeStack = null;

  protected AnatKey _key = null;
  protected AnatomyElement _anatomyArr[];
  protected keyToControllerAdaptor K2C_1 = null;
  protected int _nAnatRows;

  private HelpSet _hs = null;
  private HelpBroker _helpBroker = null;
  private final boolean _debug = false;

  public boolean _OBJModelReady = false;
  public boolean _anatomyBuilt = false;
  public boolean _greyLevelClosed = false;

  private String SLASH = System.getProperty("file.separator");

  protected static int defViewW = 450;
  protected static int defViewH = 550;

//-----------------------------------------------------
  /**
   * Purpose:    Constructor
   **/
  public SVParent2D() {
    _Locks = new SVLocks();
    init2D();
    initHelp();
  }

  public SVParent2D(File wlzFile) {
    super();
    openGreyLevel(wlzFile);

  }

//-----------------------------------------------------
  public SVLocks getSVLocks() {
     return _Locks;
  }
//-----------------------------------------------------
  private final void root() {
     /* only exists so that reflection can determine
        if this is SVParent2D or a sub-class */
  }
//-----------------------------------------------------
  protected void init2D() {
    _openViews = new Vector();
    _key = AnatKey.instance();
    _key.setSize(900,500);
    _key.pack();
    _key.setVisible(false);
    _key.setResizable(true);
    _nAnatRows = AnatKey._nrows;
    _anatomyArr = new AnatomyElement[_nAnatRows];

    for (int i=0; i<_nAnatRows; i++) {
      _anatomyArr[i] = null;
    }

    if (K2C_1 != null) _key.removeActionListener(K2C_1);
    K2C_1 = new keyToControllerAdaptor();
    _key.addActionListener(K2C_1);
  }
//-----------------------------------------------------
  protected void initHelp() {

    URL hsURL = null;
    Dummy dummy = new Dummy();

    try {
      hsURL = HelpUtils.getHelpSetURL(dummy,"help/ViewerHelp.hs");
      _hs = new HelpSet(null, hsURL);
      _helpBroker = _hs.createHelpBroker();
      _helpBroker.initPresentation();
      _helpBroker.setDisplayed(false);
    } catch(Exception ex) {
      System.out.println(ex.getMessage());
      System.out.println("error initialising help");
    }
  }

//-----------------------------------------------------
  public SVFrame createExternalView(String viewstr) {

     SVFrame ret = null;
     SectionViewer sv = null;

     ret = new SVFrame(viewstr);
     sv = new SectionViewer(viewstr, this);

     if((ret != null) && (sv != null)) {
        ret.setSectionViewer(sv);
	sv.setFrame(ret);
        addView(sv);
     }
     return ret;
  }
//-----------------------------------------------------
  public SVFrameInt createInternalView(String viewstr) {

     SVFrameInt ret = null;
     SectionViewer sv = null;

     ret = new SVFrameInt(viewstr);
     sv = new SectionViewer(viewstr, this);

     if((ret != null) && (sv != null)) {
        ret.setSectionViewer(sv);
	sv.setFrame(ret);
        addView(sv);
     }
     return ret;
  }
//-----------------------------------------------------
  public void reset() {
    closeGreyLevel();
  }

//======================================================
/**
 * Purpose:
 * @param:     File
 **/
  public void openGreyLevel(File imgFile) {
  
     closeGreyLevel();

// check that it's a Wlz file
    if(!imgFile.getName().endsWith(".wlz")) return;

    synchronized(_Locks._svpLock3) {
       while(!_greyLevelClosed) {
	  try {
	     _Locks._svpLock3.wait();
	  }
	  catch (InterruptedException iex1) {
	  }
       }
    }

    try {
      synchronized(_Locks._svpLock2) {
	 _OBJModel = new WlzObjModel(imgFile);
	 if(_OBJModel != null) {
	    _OBJModelReady = true;
            _greyLevelClosed = false;
	    _Locks._svpLock2.notifyAll();
	 } else {
	    System.out.println("_OBJModel == null");
	 }
      } // sync

      StringBuffer filstr = new StringBuffer(imgFile.getAbsolutePath());
      int fullLen = filstr.length();
      int namLen = imgFile.getName().length();
      String stagestr = filstr.substring(0, fullLen-namLen-1);
      _titleText = stagestr;
      // for multithreading
      Thread t = new Thread(rBuildAnatomy);

      t.start();
    } catch(Exception e) {
      System.out.println(e.getMessage());
    }
  } // openGreyLevel()

//-----------------------------------------------------
  public void closeGreyLevel() {
    if(_debug) System.out.println("closing grey level");
    synchronized(_Locks._svpLock3) {
       _greyLevelClosed = false;
       _titleText = "no grey level file";
       _OBJModel = null;
       _OBJModelReady = false;
       _anatBuilder = null;
       _anatomyBuilt = false;
       clearAnatomy();
       _key.setVisible(false);

       SectionViewer SV = null;
       while(_openViews.isEmpty() == false) {
	 SV = (SectionViewer)_openViews.elementAt(0);
	 SV.close();
       }
       _openViews.clear();
       _greyLevelClosed = true;
       if(_debug) System.out.println("_Locks._svpLock3 notifying");
       _Locks._svpLock3.notifyAll();
    } // sync
  } // closeGreyLevel()

//-----------------------------------------------------
  public void addView(SectionViewer view) {
    _openViews.add(view);
  }

//-----------------------------------------------------
  public void removeView(SectionViewer view) {
    _openViews.remove(view);
  }
//-----------------------------------------------------
/**
 * Purpose:    getOBJModel
 * @param:     void
 **/
  public WlzObjModel getOBJModel() {
    synchronized(_Locks._svpLock2) {
       while(!_OBJModelReady) {
	  try {
	     if(_debug) System.out.println("SVParent2D waiting for OBJModel");
	     _Locks._svpLock2.wait();
	  }
	  catch (InterruptedException iex1) {
	  }
       }
       if(_debug) System.out.println("SVParent2D returning OBJModel");
       return _OBJModel;
    } // sync
  }

//-----------------------------------------------------
/**
 * Purpose:    getAnatomyBuilder
 * @param:     void
 **/
  public AnatomyBuilder getAnatomyBuilder() {
    return _anatBuilder;
  }

//-----------------------------------------------------
/**
 * Purpose:    getSVTitleText
 * @param:     void
 **/
  public String getSVTitleText() {
    return _titleText;
  }

//-----------------------------------------------------
  public Vector collectWlzFiles(String str, Vector sofar) {
    File currDir = new File(str);
    //Vector cumTotal = new Vector(sofar);
    Vector cumTotal = sofar;
    if(!currDir.isDirectory()) {
      //System.out.println("not a directory");
      return cumTotal;
    }

    File entries[] = null;
    entries = currDir.listFiles();
    int len = entries.length;

    for(int i=0; i<len; i++) {
      //System.out.println(entries[i].getName());
      if(entries[i].isDirectory()) {
        collectWlzFiles(entries[i].getAbsolutePath(), cumTotal);
      } else {
        if(isValidName(currDir.getName(), entries[i].getName())) {
          cumTotal.add(entries[i].getAbsolutePath());
        }
      }
    }
    return cumTotal;
  }

//-------------------------------------------------------------
  public void showAnatKey() {
    if(_key.isShowing()) return;
    _key.setVisible(true);
  }

//-------------------------------------------------------------
/**
 * Purpose:    clearAnatomy
 * @param:     void
 **/
  public void clearAnatomy() {
    SectionViewer SV = null;
    int numViews = _openViews.size();
    if(numViews <= 0) return;

    for(int i=0; i<numViews; i++) {
      SV = (SectionViewer)_openViews.elementAt(i);
      SV.removeAnatomy();
      SV.setAnatomyText(SV._anatomyStr);
    }

    for(int i=0; i<_nAnatRows; i++) {
      _anatomyArr[i] = null;
    }
    AnatKey.reset();
  }

//-------------------------------------------------------------
  public void showAnatomy(String str) {

    WlzFileInputStream in;
    File thisFile = new File(str);

    SectionViewer SV = null;

    int pathlen = _anatBuilder._pathLengthToAnatomy;
    String fullName = thisFile.getAbsolutePath();
    String anatName = fullName.substring(pathlen);
    int len = anatName.length() - 4; // chop off '.wlz'
    String anatName1 = anatName.substring(0,len);

    try {
      in = new WlzFileInputStream(str);
      updateAnatomyArr(WlzObject.WlzReadObj(in), capitalise(anatName1)+" *");
    } catch(WlzException we) {
      System.out.println("showAnatomy: 1");
      System.out.println(we.getMessage());
      if(we.getMessage().indexOf("WLZ_ERR_OBJECT_NULL") != -1) {
        updateAnatomyArr(null, capitalise(anatName1));
      }
    } catch (IOException ie) {
      System.out.println(ie.getMessage());
    }

    showAnatKey();

    int numViews = _openViews.size();
    if(numViews == 0) return;
    for(int i=0; i<numViews; i++) {
      SV = (SectionViewer)_openViews.elementAt(i);
      SV.removeAnatomy();
      SV.anatomyFromMenu(_anatomyArr);
    }
  }

//-------------------------------------------------------------
  public void showCombinedAnatomy(String str) {
//    Vector files = null; // String pathname to file or dir

    File thisDir = new File(str);
    Vector files = collectWlzFiles(str, new Vector());

    int pathlen = _anatBuilder._pathLengthToAnatomy;
    String fullName = thisDir.getAbsolutePath();
    String anatName = fullName.substring(pathlen);
  
    //System.out.println("anatName = "+anatName);

    try {
      updateAnatomyArr(combineWlzObjs(files), capitalise(anatName));
    } catch (Exception e) {
      System.out.println("Can not load the anatomy!");
    }

    showAnatKey();

    SectionViewer SV = null;

    int numViews = _openViews.size();
    if(numViews <= 0) return;

    for(int i=0; i<numViews; i++) {
      SV = (SectionViewer)_openViews.elementAt(i);
      SV.removeAnatomy();
      SV.anatomyFromMenu(_anatomyArr);
    }
  }

//-------------------------------------------------------------
/**
 * Purpose:    combineWlzObjs
 * @param:     Vector
 **/
  public WlzObject combineWlzObjs(Vector files) {
    WlzObject unionObj = null;
    WlzObject obj1 = null;
    WlzObject obj2 = null;
    WlzFileInputStream in = null;

    int len = files.size();
    if(len == 0) return (WlzObject)null;

    try {
      in = new WlzFileInputStream((String)files.elementAt(0));
      obj1 = WlzObject.WlzReadObj(in);
    } catch(WlzException we) {
      System.out.println("combineWlzObjs: 1");
      System.out.println(we.getMessage());
      if(we.getMessage().indexOf("WLZ_ERR_OBJECT_NULL") != -1) {
        obj1 = null;
      }
    } catch (IOException ie) {
      System.out.println(ie.getMessage());
    }

    if(len == 1) return obj1;

    if(obj1 != null) {
      try {
        for(int i=1; i<len; i++) {
          in = new WlzFileInputStream((String)files.elementAt(i));
          obj2 = WlzObject.WlzReadObj(in);
          if (obj2 != null) unionObj = WlzObject.WlzUnion2(obj1, obj2);

          in = null;
          obj1 = null;
          obj2 = null;
          obj1 = unionObj;
        }
      } catch(WlzException we) {
        System.out.println("combineWlzObjs: 2");
        System.out.println(we.getMessage());
      } catch (IOException ie) {
        System.out.println(ie.getMessage());
      }
    }
    return unionObj;
  }

//-------------------------------------------------------------
  public void updateAnatomyArr(WlzObject obj, String str) {
   /*
      the array elements will be filled
      as items of anatomy are selected.
      If more anatomy items are selected
      than the length of the array,
      the first one(s) are overwritten

      a new item occupies the lowest available position
    */

    String descr = str;
    if(obj == null) descr = new String(str+" (not painted)");

    int indx = AnatomyElement.getNextIndex(_anatomyArr);
    _anatomyArr[indx] = new AnatomyElement(obj, descr, indx);
    AnatKey.update(_anatomyArr);
  }

//-------------------------------------------------------------
  protected boolean isValidName(String dir, String fil) {
    if(!fil.endsWith(".wlz")) return false;
    int len = fil.length() - 4; // chop off '.wlz'
    if((fil.substring(0,len)).equals(dir)) {
      return true;
    } else {
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
    if(lastSlash != -1) {
       endBit = buf.substring(lastSlash);
       buf.replace(lastSlash, buf.length(), endBit.toUpperCase());
       ret = buf.toString();
    } else if(name.equals("embryo") ||
              name.equals("extraembryonic_component")) {
       ret = name.toUpperCase();
    }

    return ret;
  }

//-------------------------------------------------------------
  public AnatKey getAnatomyKey() {
    return _key;
  }

//-------------------------------------------------------------
  public AnatomyElement[] getAnatomyArr() {
    return _anatomyArr;
  }

//-------------------------------------------------------------
  public HelpSet getHelpSet() {
    return _hs;
  }

//-------------------------------------------------------------
  public HelpBroker getSVHelpBroker() {
    return _helpBroker;
  }

//-------------------------------------------------------------
  public Vector getOpenViews() {
    return _openViews;
  }

//===============================================================
// Thread stuff
//===============================================================

  Runnable rBuildAnatomy = new Runnable() {
     public void run() {

	try {
	   _anatBuilder = new AnatomyBuilder();
	   _anatBuilder.buildAnatomy(new File(_titleText));
	   synchronized(_Locks._svpLock1) {
	      _anatomyBuilt = true;
	      if(_debug) System.out.println("_Locks._svpLock1 notifying");
	      _Locks._svpLock1.notifyAll();
	   } // sync
	}
	catch(Exception e) {
	   System.out.println("SVParent2D rBuildAnatomy thread: exception");
	   e.printStackTrace();
	}
     }
  };



//===============================================================
  public class keyToControllerAdaptor implements ActionListener {

    int numViews;
    int indx;
    SectionViewer SV = null;
    Color newCol = null;

    public keyToControllerAdaptor() {
    }

    public void actionPerformed(ActionEvent e) {
      String cmd = e.getActionCommand();
      boolean wasRemoved = false;

      indx = (
          Integer.valueOf(cmd.substring(cmd.length()-1))).intValue();
      //if(_anatomyArr[indx] == null) return;

      if(cmd.indexOf("makeVisible") != -1) {
	if(_anatomyArr[indx] != null) {
	   _anatomyArr[indx].setVisible(true);
	}
      } else if(cmd.indexOf("makeInvisible") != -1) {
	if(_anatomyArr[indx] != null) {
	   _anatomyArr[indx].setVisible(false);
	}
      } else if(cmd.indexOf("zap") != -1) {
	if(_anatomyArr[indx] != null) {
	   wasRemoved = _anatomyArr[indx].isRemoved();
	   _key.setZapToolTips(indx, !wasRemoved);
	   _anatomyArr[indx].setRemoved(!wasRemoved);
	}
      } else if(cmd.indexOf("colour") != -1) {
         newCol = JColorChooser.showDialog(null,
	                                   "choose colour for anatomy component",
					   new Color(_key._cols[indx]));
         _key.setCol(newCol, indx);
      } else {
        System.out.println("unknown command");
      }

      AnatKey.update(_anatomyArr);
      doChange();
    }

    private void doChange() {
      numViews = _openViews.size();
      if(numViews == 0) return;

      for(int i=0; i<numViews; i++) {
        SV = (SectionViewer)_openViews.elementAt(i);
        SV.removeAnatomy();
        SV.anatomyFromMenu(_anatomyArr);
      }
    }
  }

//-------------------------------------------------------------
}
