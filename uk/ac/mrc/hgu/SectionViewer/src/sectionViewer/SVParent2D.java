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

/**
 *   Utility class which manages SectionViewers for an application.
 *   If an application doesn't use SVParent2D then one of its classes
 *   <b>must</b> implement SVParent.
 */
public class SVParent2D implements SVParent {

  /**   Directory containing the underlying 3D Woolz object */
  private String _wlzDir = "";

  /**   Model for the underlying 3D Woolz object. */
  private WlzObjModel _OBJModel = null;
  
  /**   Object that will build the anatomy menus. */
  private AnatomyBuilder _anatBuilder = null;

  /**   Object that provides locks for multi-threading synchronization. */
  private SVLocks _Locks = null;

  /**   Collection of currently open SectionViewers. */
  protected Vector _openViews = null;

  /**   Indicates colour of anatomy components in SectionViewers. */
  protected AnatKey _key = null;

  /**
   * Array of objects corresponding to rows in the AnatKey.
   * Obsolete.
   */
  protected AnatomyElement _anatomyArr[];

  /**   Array of objects corresponding to rows in the AnatKey. */
  protected Vector _elementVec = null;

  /**   Manages events initiated from the AnatKey controls. */
  protected keyToControllerAdaptor K2C_1 = null;


  /**   Java Help Helpset for SectionViewer. */
  private HelpSet _hs = null;

  /**   Java Help HelpBroker for SectionViewer. */
  private HelpBroker _helpBroker = null;

  /**   toggles display of debugging messages */
  private final boolean _debug = false;

  /**   Indicates that the 3D Woolz object has been read in */
  public boolean _OBJModelReady = false;

  /**   Indicates that the anatomy menus have been built */
  public boolean _anatomyBuilt = false;

  /**   Indicates that the Woolz object has been closed */
  public boolean _greyLevelClosed = false;

  /**   The system file separator ('\' or '/') */
  private String SLASH = System.getProperty("file.separator");

  /**   Default width of SectionViewer */
  protected static int defViewW = 450;

  /**   Default height of SectionViewer */
  protected static int defViewH = 550;

//-----------------------------------------------------
  /**   Creates an SVParent2D. */
  public SVParent2D() {
    _Locks = new SVLocks();
    init2D();
    initHelp();
  }

  /**
   *   Creates an SVParent2D and opens the given 3D Woolz file.
   *   @param wlzFile the 3D Woolz file to open.
   */
  public SVParent2D(File wlzFile) {
    super();
    openGreyLevel(wlzFile);

  }

//-----------------------------------------------------
  /**
   *   Returns the SVLocks object containing locks required for synchronization.
   *   @return the instance of SVLocks.
   */
  public SVLocks getSVLocks() {
     return _Locks;
  }
//-----------------------------------------------------
  /**
   *   Exists so that SectionViewers can determine by Reflection
   *   if this is an SVParent2D or a sub-class
   */
  private final void root() {
  }
//-----------------------------------------------------
  /**
   *   Initialises the collection of open SectionViewers
   *   and the AnatKey.
   */
  protected void init2D() {
    _openViews = new Vector();
    /* make sure only 1 AnatKey is created */
    if(_key == null) {
       _key = new AnatKey();
       _key.setTitle("Anatomy Key");
    }
    _key.pack();
    _key.setVisible(false);
    _key.setResizable(true);
    _elementVec = new Vector();

    K2C_1 = new keyToControllerAdaptor();
  }
//-----------------------------------------------------
  /**
   *   Initialises the Java Help system for SectionViewer.
   */
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
  /**
   *   Tests for the existence of an AnatKey.
   *   @return true if there is an AnatKey.
   */
  public boolean hasAnatKey() {
     return (_key == null) ? false : true;
  }
//-----------------------------------------------------
  /**
   *   Factory method for SectionViewers in external windows.
   *   @param viewstr the initial orientation of the SectionViewer.
   *   @return SVFrame containing a new SectionViewer of the given orientation.
   */
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
  /**
   *   Factory method for SectionViewers in internal windows.
   *   @param viewstr the initial orientation of the SectionViewer.
   *   @return SVFrameInt containing a new SectionViewer of the given orientation.
   */
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
  /**
   *   Closes the underlying 3D Woolz object.
   */
  public void reset() {
    closeGreyLevel();
  }

//======================================================
  /**
   *   Opens a 3D Woolz object.
   *   @param imgFile the 3D Woolz file.
   */
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
      _wlzDir = stagestr;
      // for multithreading
      Thread t = new Thread(rBuildAnatomy);

      t.start();
    } catch(Exception e) {
      System.out.println(e.getMessage());
    }
  } // openGreyLevel()

//-----------------------------------------------------
  /**
   *   Closes the underlying 3D Woolz object.
   */
  public void closeGreyLevel() {
    if(_debug) System.out.println("closing grey level");
    synchronized(_Locks._svpLock3) {
       _greyLevelClosed = false;
       //_wlzDir = "no grey level file";
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
  /**
   *   Adds the given SectionViewer to the collection
   *   of open SectionViewers.
   *   @param view the SectionViewer to add.
   */
  public void addView(SectionViewer view) {
    _openViews.add(view);
  }

//-----------------------------------------------------
  /**
   *   Removes the given SectionViewer from the collection
   *   of open SectionViewers.
   *   @param view the SectionViewer to remove.
   */
  public void removeView(SectionViewer view) {
    _openViews.remove(view);
  }
//-----------------------------------------------------
  /**
   *   Returns the Model for the underlying 3D Woolz object.
   *   This is synchronized so that it can't return until 
   *   the 3D Woolz object has been read in.
   *   @return the WlzObjModel for the current 3D Woolz object.
   */
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
   *   Returns an object used to build anatomy menus.
   *   @return a new AnatomyBuilder.
   */
  public AnatomyBuilder getAnatomyBuilder() {
    return _anatBuilder;
  }

//-----------------------------------------------------
  /**
   *   Obsolete
   */
  public String getSVTitleText() {
    //return _wlzDir;
    return "";
  }

//-----------------------------------------------------
  /**
   *   Walks through an anatomy directory hierarchy accumulating 
   *   3D Woolz files for each available anatomy component.
   *   @param str the filename of the current directory.
   *   @param sofar the current collection of anatomy component
   *   pathnames.
   *   @return the collection of anatomy component pathnames.
   */
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
  /**
   *   Makes the AnatKey visible.
   */
  public void showAnatKey() {
    if(_key.isShowing()) return;
    _key.setVisible(true);
  }

//-------------------------------------------------------------
  /**
   *   Returns the AnatKey to its initial empty state.
   */
  public void clearAnatomy() {
    SectionViewer SV = null;
    int numViews = _openViews.size();
    if(numViews <= 0) return;

    for(int i=0; i<numViews; i++) {
      SV = (SectionViewer)_openViews.elementAt(i);
      SV.removeAnatomy();
      SV.setAnatomyText(SV._anatomyStr);
    }

    _elementVec.clear();
    _key.reset();
  }

//-------------------------------------------------------------
  /**
   *   Adds the given (atomic) component to the collection
   *   of anatomy components and causes it to be displayed
   *   in all open SectionViewers.
   *   @param str the filename of the anatomy component to add.
   */
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
      updateElementVec(WlzObject.WlzReadObj(in),
                       capitalise(anatName1)+" *");
    }
    catch(WlzException we) {
      System.out.println("showAnatomy: 1");
      System.out.println(we.getMessage());
      if(we.getMessage().indexOf("WLZ_ERR_OBJECT_NULL") != -1) {
        updateElementVec(null, capitalise(anatName1));
      }
    }
    catch (IOException ie) {
      System.out.println(ie.getMessage());
    }

    showAnatKey();

    int numViews = _openViews.size();
    if(numViews == 0) return;
    for(int i=0; i<numViews; i++) {
      SV = (SectionViewer)_openViews.elementAt(i);
      SV.removeAnatomy();
      SV.anatomyFromMenu(_elementVec);
    }
  }

//-------------------------------------------------------------
  /**
   *   Adds the given (high level) component to the collection
   *   of anatomy components and causes it to be displayed
   *   in all open SectionViewers.
   *   @param str the filename of the anatomy component to add.
   */
  public void showCombinedAnatomy(String str) {

    File thisDir = new File(str);
    Vector files = collectWlzFiles(str, new Vector());

    int pathlen = _anatBuilder._pathLengthToAnatomy;
    String fullName = thisDir.getAbsolutePath();
    String anatName = fullName.substring(pathlen);
  
    try {
      updateElementVec(combineWlzObjs(files),
                       capitalise(anatName));
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
      SV.anatomyFromMenu(_elementVec);
    }
  }

//-------------------------------------------------------------
  /**
   *   Combines atomic components into one high level component.
   *   @param files the collection of atomic component filenames.
   *   @return the high level WlzObject.
   */
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
  /**
   *   Updates the Vector of AnatomyElements represented by
   *   the rows in the AnatKey.
   *   @param obj 3D Woolz object representing an anatomy component.
   *   @param str the full path name of the anatomy component.
   */
  public void updateElementVec(WlzObject obj, String str) {

    int indx = 0;

    String descr = str;
    if(obj == null) descr = new String(str+" (not painted)");

    indx = _key.nextIndx();
    AnatomyElement el = new AnatomyElement(obj, descr, indx);
    _elementVec.add(el);

    _key.addRow(descr);
    KeyEntry row = _key.getRow(indx);
    if(row != null) {
       row.addActionListener(K2C_1);
    }
  }

//-------------------------------------------------------------
  /**
   *   Checks that the given filename is valid.
   *   @param dir the directory containing the given file.
   *   @param fil the full path name of the given file.
   *   @return true if the given filename ends in '.wlz'
   *   and the base of the name is the same as the 
   *   directory that contains it.
   */
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
  /**
   *   Capitalises the final part of the given filename.
   *   @param name a full filename.
   *   @return the given filename with the final part capitalised.
   */
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
  /**
   *   Returns the AnatKey.
   *   @return _key.
   */
  public AnatKey getAnatomyKey() {
    return _key;
  }

//-------------------------------------------------------------
  /**
   *   Returns the collection of objects represented by
   *   the rows of the AnatKey.
   *   @return the AnatomyElement Vector.
   */
  public Vector getAnatomyElements() {
    return _elementVec;
  }

//-------------------------------------------------------------
  /**
   *   Returns the Java Help Helpset for SectionViewer.
   *   @return the Java Help Helpset for SectionViewer.
   */
  public HelpSet getHelpSet() {
    return _hs;
  }

//-------------------------------------------------------------
  /**
   *   Returns the Java Help Help broker for SectionViewer.
   *   @return the Java Help Help broker for SectionViewer.
   */
  public HelpBroker getSVHelpBroker() {
    return _helpBroker;
  }

//-------------------------------------------------------------
  /**
   *   Returns the collection of currently open SectionViewers.
   *   @return the _openViews Vector.
   */
  public Vector getOpenViews() {
    return _openViews;
  }

//-------------------------------------------------------------
  /**
   *   Returns the AnatomyElement with the specified indx.
   *   @param indx a unique index.
   *   @return the AnatomyElement with the specified indx.
   */
  public AnatomyElement getAnatomyElement(int indx) {

     Enumeration els = _elementVec.elements();
     AnatomyElement el = null;
     while (els.hasMoreElements()) {
        el = (AnatomyElement)els.nextElement();
        if(el.getIndx() == indx) break;
     }
     return el;
  }


//===============================================================
// Thread stuff
//===============================================================

  /**   Builds the anatomy menus in a new thread.**/
  Runnable rBuildAnatomy = new Runnable() {
     public void run() {

	try {
	   _anatBuilder = new AnatomyBuilder();
	   _anatBuilder.buildAnatomy(new File(_wlzDir));
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
  /**
   *   Listens for ActionEvents from an AnatKey and implements the
   *   appropriate action depending upon the ActionCommand.
   */
  public class keyToControllerAdaptor
      implements ActionListener {

    int numViews;
    int indx;
    SectionViewer SV = null;

    public keyToControllerAdaptor() {
    }

    public void actionPerformed(ActionEvent e) {
      String cmd = e.getActionCommand();

      String[] splits = cmd.split("\\d", 2);
      /* 
       * assumes cmd is of the form "name123"
       * where 'name' describes the action to be taken 
       * and '123' is the unique index number of the keyEntry
       */
      indx = (Integer.valueOf(
                cmd.substring(splits[0].length()))).intValue();

      if(cmd.indexOf("makeVisible") != -1) {
         _key.setEntryVisible(indx, true);
	 getAnatomyElement(indx).setElementVisible(true);
      } else if(cmd.indexOf("makeInvisible") != -1) {
         _key.setEntryVisible(indx, false);
	 getAnatomyElement(indx).setElementVisible(false);
      } else if(cmd.indexOf("zap") != -1) {
	 _key.removeRow(indx);
	 _elementVec.remove(getAnatomyElement(indx));
      } else if(cmd.indexOf("colour") != -1) {
         _key.setCol(indx);
      } else {
        System.out.println("unknown command");
      }

      doChange();
    }

    private void doChange() {
      numViews = _openViews.size();
      if(numViews == 0) return;

      for(int i=0; i<numViews; i++) {
        SV = (SectionViewer)_openViews.elementAt(i);
        SV.removeAnatomy();
        SV.anatomyFromMenu(_elementVec);
      }
    }
  }

//-------------------------------------------------------------
}
