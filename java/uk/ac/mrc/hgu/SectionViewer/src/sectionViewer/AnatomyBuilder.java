package sectionViewer;
import sectionViewer.*;

import javax.swing.*;
import java.util.*;
import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

/**
 *   Builds anatomy component menus by inspecting
 *   the anatomy directory hierarchy.
 *   Multi-threading is used.
 */
public class AnatomyBuilder {

  /**
   *   Anatomy component file names in the
   *   'embryo' hierarchy
   */
  private volatile Stack _embFileStack = null;

  /**
   *   Anatomy component file names in the
   *   'extra_embryonic_component' hierarchy
   */
  private volatile Stack _xembFileStack = null;
  
  /**
   *   Anatomy component Woolz objects in the
   *   'embryo' hierarchy
   */
  private volatile Stack _embObjStack = null;

  /**
   *   Anatomy component Woolz objects in the
   *   'extra_embryonic_component' hierarchy
   */
  private volatile Stack _xembObjStack = null;

  /**
   *   Flag to indicate when anatomy menus have been built.
   *   Tested by multi-threading code.
   */
  private boolean _done;

  /**
   *   System file separator ('/' or '\')
   */
  private String SLASH = System.getProperty("file.separator");

  /**
   *   Menu for anatomy components from the 'embryo' hierarchy.
   */
  SelectableMenu _embMenu = null;

  /**
   *   Menu for anatomy components from
   *   the 'extra_embryonic_component' hierarchy.
   */
  SelectableMenu _xembMenu = null;

  /**
   *   The path length up to the final part of the anatomy component name.
   */
  protected static int _pathLengthToAnatomy = 0;

  /**
   *   Prevents access to anatomy menus until they have been built.
   */
  public final Object _lock = new Object();

  /**
   *   Returns the lock object which prevents access to anatomy menus
   *   until they have been built.
   *   @return _lock 
   */
  public Object getLock() {
     return _lock;
  }

  // constructor
  /**
   *   Default constructor
   */
  public AnatomyBuilder() {
     _embFileStack = new Stack();
     _xembFileStack = new Stack();
     _embObjStack = new Stack();
     _xembObjStack = new Stack();
     _embMenu = new SelectableMenu("embryo");
     _xembMenu = new SelectableMenu("extra_embryonic_component");
     _done = false;
  }

//-------------------------------------------------------------
// methods for walking the anatomy directory
//-------------------------------------------------------------
  /**
   *   Organises building of anatomy menus.
   *   @param stageDir is the full name of the top-level directory for 
   *   a structure containing anatomy components. For example with mouse embryos
   *   this would be the top level directory for a Theiler stage, such as ts14.
   */
  public void buildAnatomy(File stageDir) {
     _done = false;
     File anaDir = null;
     File embDir = null;
     File xembDir = null;
     Stack fileStack = new Stack();

     if(!stageDir.isAbsolute()) return;

     synchronized(_lock) {
	anaDir = new File(stageDir.getAbsolutePath() + SLASH + "anatomy");
	if(anaDir.exists() != true) return;

	_pathLengthToAnatomy = anaDir.getAbsolutePath().length() + 1;

	embDir = new File(anaDir.getAbsolutePath() + SLASH + "embryo");
	xembDir = new File(anaDir.getAbsolutePath() + SLASH + "extraembryonic_component");

	collectWlzFiles(embDir, _embFileStack);
	collectWlzFiles(xembDir, _xembFileStack);
	makeObjStack(0); // embryonic
	makeObjStack(1); // extra-embryonic
	buildMenu(embDir, _embMenu);
	buildMenu(xembDir, _xembMenu);
	_done = true;
	//printAnatomy();
	_lock.notifyAll();
     } // sync
  }

//----------------------------------------------------
  /**
   *   Recursive function which collects Woolz files.
   *   @param dir the directory in which to look for Woolz files.
   *   If the directory contains sub-directories, these are searched as well.
   *   @param stack the Stack on which to store Woolz filenames.
   */
  protected void collectWlzFiles(File dir, Stack stack) {

     // walk through the directories
     String fullPath = dir.getAbsolutePath();
     String contents[] = null;
     File thisFile = null;

     if(!dir.isDirectory()) return;

     contents = dir.list();

     int len = contents.length;

     for (int i=0; i < len; i++) {
	thisFile = new File(fullPath + SLASH + contents[i]);
	if (thisFile.isDirectory()) {
	   //System.out.println("dir = "+thisFile.getName());
	   collectWlzFiles(thisFile, stack);
	} else if (thisFile.isFile()) {
	   if (thisFile.getName().lastIndexOf(".wlz") > 0) {
	      stack.push(thisFile);
	   }
	}
     }
  } // collectWlzFiles()

//----------------------------------------------------
  /**
   *   Makes a Stack of Woolz Objects from anatomy component filenames.
   *   @param type = 0 for embryonic, 1 for extra_embryonic_component.
   */
  protected void makeObjStack(int type) {

     File thisFile = null;
     WlzObject obj = null;
     WlzFileInputStream in = null;

     Stack fstack = null;
     Stack ostack = null;

     int len = 0;

     switch(type){
        case 0: // embryonic
	   fstack = _embFileStack;
	   ostack = _embObjStack;
	   break;
        case 1: // extra-embryonic
	   fstack = _xembFileStack;
	   ostack = _xembObjStack;
	   break;
        default:
	   break;
     }

     len = fstack.size();

     try {
	for (int i=0; i<len; i++) {
	   thisFile = (File)fstack.elementAt(i);
	   in = new WlzFileInputStream(thisFile.getAbsolutePath());
	   obj = WlzObject.WlzReadObj(in);
	   ostack.push(obj);
	}
     }
     catch (WlzException e) {
        System.out.println("makeObjStack");
	System.out.println(e.getMessage());
     }
     catch (IOException e) {
	System.out.println(e.getMessage());
     }

  } // makeObjStack()

//----------------------------------------------------
  /**
   *   Recursive function which builds a SelectableMenu.
   *   @param dir the directory in which to look for Woolz file names.
   *   If the directory contains sub-directories, these are searched as well.
   *   @param menu the Selectable menu.
   */
  protected void buildMenu(File dir, SelectableMenu menu) {
     // walk through the directories
     String fullPath = dir.getAbsolutePath();
     menu.setActionCommand(fullPath);
     String contents[] = null;
     File thisFile = null;

     if(!dir.isDirectory()) return;

     contents = dir.list();

     int len = contents.length;

     for (int i=0; i<len; i++) {
	thisFile = new File(fullPath + SLASH + contents[i]);
	if (thisFile.isDirectory()) {
	   SelectableMenu subMenu = new SelectableMenu(thisFile.getName());
	   subMenu.setActionCommand(fullPath + SLASH + contents[i]);
	   menu.add(subMenu);

	   int newlen = thisFile.list().length;
	   if(newlen > 0) {
	      buildMenu(thisFile, subMenu);
	   }
	} else {
	   if(thisFile.getName().equals(
	      thisFile.getParentFile().getName() + ".wlz")) {

	      JMenuItem item = new JMenuItem(
	               leafMod(thisFile.getName()));
	      item.setActionCommand(fullPath + SLASH + contents[i]);
	      menu.add(item, 0);
	   }
	}
     }
  }

//----------------------------------------------------
  /**
   *   Attaches the string "(domain)" to menu entries
   *   if they are atomic components.
   *   @param nam the name of the atomic component.
   */
  private String leafMod(String nam) {

     StringBuffer buf = new StringBuffer(nam);
     int len = buf.length();
     return (buf.replace(len-4, len, " (domain)")).toString();
  }

//----------------------------------------------------
  /**
   *    Returns true when anatomy component menus have been built.
   *    @return _done
   */
  public boolean getDone() {
     return _done;
  }

//----------------------------------------------------
  /**
   *    Utility test function
   */
  public void printAnatomy() {
     Enumeration lmnts = null;
      try {
         lmnts = _embFileStack.elements();
         System.out.println("EMBRYO");
         while(lmnts.hasMoreElements() == true) {
            System.out.println(((File)lmnts.nextElement()).getName());
         }
         lmnts = _xembFileStack.elements();
         System.out.println("XEMBRYO");
         while(lmnts.hasMoreElements() == true) {
            System.out.println(((File)lmnts.nextElement()).getName());
         }
      }
      catch(NoSuchElementException e) {
         System.out.println(e.getMessage());
      }
  }

//----------------------------------------------------
  public Stack getEmbFileStack() { return _embFileStack; }
  public Stack getXEmbFileStack() { return _xembFileStack; }
  public Stack getEmbObjStack() { return _embObjStack; }
  public Stack getXEmbObjStack() { return _xembObjStack; }

  public SelectableMenu getEmbMenu() { return _embMenu; }
  public SelectableMenu getXEmbMenu() { return _xembMenu; }

  public int getPathLengthToAnatomy() { return _pathLengthToAnatomy; }

} // class AnatomyBuilder
