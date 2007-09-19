package sectionViewer;

import java.io.*;
import java.util.*;
import java.util.List;

import java.awt.*;
import javax.swing.*;

import uk.ac.mrc.hgu.Wlz.*;

/**
 *   Builds anatomy component menus by inspecting
 *   the anatomy directory hierarchy.
 *   Multi-threading is used.
 */
public class AnatomyBuilder {

  /**
   *   Anatomy component file names
   */
  private volatile Stack<File> _anaFileStack = null;

  /**
   *   Anatomy component file names in the
   *   'embryo' hierarchy
   */
  private volatile Stack<File> _embFileStack = null;

  /**
   *   Anatomy component file names in the
   *   'extra_embryonic_component' hierarchy
   */
  private volatile Stack<File> _xembFileStack = null;

  /**
   *   Anatomy component Woolz objects
   */
  private volatile Stack<WlzObject> _anaObjStack = null;

  /**
   *   Anatomy component Woolz objects in the
   *   'embryo' hierarchy
   */
  private volatile Stack<WlzObject> _embObjStack = null;

  /**
   *   Anatomy component Woolz objects in the
   *   'extra_embryonic_component' hierarchy
   */
  private volatile Stack<WlzObject> _xembObjStack = null;

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
   *   Menu for anatomy components.
   */
  SelectableMenu _anaMenu = null;

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

  //private Font _menuFont = new Font("default", Font.PLAIN, 11);
  private Font _menuFont = null;

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
   *   Constructs an AnatomyBuilder with the default Font.
   */
  public AnatomyBuilder() {
     this(null);
  }
  /**
   *   Constructs an AnatomyBuilder with the given Font.
   */
  public AnatomyBuilder(Font font) {
     _anaFileStack = new Stack<File>();
     _embFileStack = new Stack<File>();
     _xembFileStack = new Stack<File>();
     _anaObjStack = new Stack<WlzObject>();
     _embObjStack = new Stack<WlzObject>();
     _xembObjStack = new Stack<WlzObject>();
     _anaMenu = new SelectableMenu("Anatomy");
     if(font == null) {
	_menuFont = _anaMenu.getFont();
     } else {
        _menuFont = font;
     }
     _anaMenu.setFont(_menuFont);
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

     //File embDir = null;
     //File xembDir = null;

     Stack<File> fileStack = new Stack<File>();

     //using anatomy tree
     if(!stageDir.isAbsolute()) return;

     synchronized(_lock) {
	anaDir = new File(stageDir.getAbsolutePath() + SLASH + "anatomy");
	if(anaDir.exists() != true) return;

	_pathLengthToAnatomy = anaDir.getAbsolutePath().length() + 1;

	//embDir = new File(anaDir.getAbsolutePath() + SLASH + "embryo");
	//xembDir = new File(anaDir.getAbsolutePath() + SLASH + "extraembryonic_component");

	// operations now combined in the 1 walk through
	//collectWlzFiles(embDir, _embFileStack);
	//collectWlzFiles(xembDir, _xembFileStack);
	//makeObjStack(0); // embryonic
	//makeObjStack(1); // extra-embryonic
	//buildMenu(embDir, _embMenu);
	//buildMenu(xembDir, _xembMenu);

	buildMenu(anaDir, _anaFileStack);
//buildMenu(anaDir, _anaMenu, _anaFileStack);
	makeObjStack(); // embryonic

	//buildMenu(embDir, _embMenu, _embFileStack);
	//buildMenu(xembDir, _xembMenu, _xembFileStack);
	//makeObjStack(0); // embryonic
	//makeObjStack(1); // extra-embryonic

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
  protected void collectWlzFiles(File dir, Stack<File> stack) {

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
   */
  protected void makeObjStack() {

     File thisFile = null;
     WlzObject obj = null;
     WlzFileInputStream in = null;

     Stack<File> fstack = null;
     Stack<WlzObject> ostack = null;

     int len = 0;

     fstack = _anaFileStack;
     ostack = _anaObjStack;

     len = fstack.size();

     try {
	for (int i=0; i<len; i++) {
	   thisFile = (File)fstack.elementAt(i);
	   in = new WlzFileInputStream(thisFile.getAbsolutePath());
	   //System.out.println("calling WlzReadObj from AnatomyBuilder ^^^^^^^^^^");
	   obj = WlzObject.WlzReadObj(in);
           in.close();
	   //System.out.println("stream has been closed ^^^^^^^^^^");
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
   *   Makes a Stack of Woolz Objects from anatomy component filenames.
   *   @param type = 0 for embryonic, 1 for extra_embryonic_component.
   */
  protected void makeObjStack(int type) {

     File thisFile = null;
     WlzObject obj = null;
     WlzFileInputStream in = null;

     Stack<File> fstack = null;
     Stack<WlzObject> ostack = null;

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
	   //System.out.println("calling WlzReadObj from AnatomyBuilder ^^^^^^^^^^");
	   obj = WlzObject.WlzReadObj(in);
           in.close();
	   //System.out.println("stream has been closed ^^^^^^^^^^");
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
     //menu.setActionCommand(fullPath);
     String[] contents = null;
     File thisFile = null;

     if(!dir.isDirectory()) return;

     contents = dir.list();

     List<String> l = Arrays.asList(contents);
     Collections.sort(l);
     contents = (String[])l.toArray();

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

  protected void buildMenu(File dir) {
     // walk through the directories
     String fullPath = dir.getAbsolutePath();
     //menu.setActionCommand(fullPath);
     String contents[] = null;
     File thisFile = null;

     if(!dir.isDirectory()) return;

     contents = dir.list();

     List<String> l = Arrays.asList(contents);
     Collections.sort(l);
     contents = (String[])l.toArray();

     int len = contents.length;

     for (int i=0; i<len; i++) {
        thisFile = new File(fullPath + SLASH + contents[i]);
        if (thisFile.isDirectory()) {
           int newlen = thisFile.list().length;
           if(newlen > 0) {
              buildMenu(thisFile);
           }
        } else {
           if(thisFile.getName().equals(
              thisFile.getParentFile().getName() + ".wlz")) {
           }
        }
     }
  }

//----------------------------------------------------
  /**
   *   Recursive function which builds a SelectableMenu and
   *   collects Woolz file names by walking the anatomy directory.
   *   @param dir the directory in which to look for Woolz file names.
   *   If the directory contains sub-directories, these are searched as well.
   *   @param menu the Selectable menu.
   *   @param stack the Stack on which to store Woolz filenames.
   */
  protected void buildMenu(File dir, SelectableMenu menu, Stack<File> stack) {
     // walk through the directories
     String fullPath = dir.getAbsolutePath();
     menu.setActionCommand(fullPath);
     String contents[] = null;
     File thisFile = null;

     if(!dir.isDirectory()) return;

     contents = dir.list();

     List<String> l = Arrays.asList(contents);
     Collections.sort(l);
     contents = (String[])l.toArray();

     int len = contents.length;

     for (int i=0; i<len; i++) {
	thisFile = new File(fullPath + SLASH + contents[i]);
	if (thisFile.isDirectory()) {
	   SelectableMenu subMenu = new SelectableMenu(thisFile.getName());
           subMenu.setFont(_menuFont);
           subMenu.setActionCommand(fullPath + SLASH + contents[i]);
	   menu.add(subMenu);

	   int newlen = thisFile.list().length;
	   if(newlen > 0) {
	      buildMenu(thisFile, subMenu, stack);
	   }
	} else {
	   if(thisFile.getName().equals(
	      thisFile.getParentFile().getName() + ".wlz")) {

	      JMenuItem item = new JMenuItem(
	               leafMod(thisFile.getName()));
                item.setFont(_menuFont);
	      item.setActionCommand(fullPath + SLASH + contents[i]);
	      menu.add(item, 0);

	      stack.push(thisFile);
	   }
	}
     }
  }

  protected void buildMenu(File dir, Stack<File> stack) {
   // walk through the directories
   String fullPath = dir.getAbsolutePath();
   String contents[] = null;
   File thisFile = null;

   if(!dir.isDirectory()) return;

   contents = dir.list();

   List<String> l = Arrays.asList(contents);
   Collections.sort(l);
   contents = (String[])l.toArray();

   int len = contents.length;

   for (int i=0; i<len; i++) {
      thisFile = new File(fullPath + SLASH + contents[i]);
      if (thisFile.isDirectory()) {
         int newlen = thisFile.list().length;
         if(newlen > 0) {
            buildMenu(thisFile, stack);
         }
      } else {
         if(thisFile.getName().equals(
            thisFile.getParentFile().getName() + ".wlz")) {
            stack.push(thisFile);
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
   *    Adds the top sub-elements to a menu
   */
  public void addSubElements(SelectableMenu menu) {

     JMenuItem item = null;

     int nitems = _anaMenu.getItemCount();

     System.out.println("# of menu items = "+nitems);

     for(int i=0; i<nitems; i++) {
        item = _anaMenu.getItem(i);
	System.out.println("menu item = "+item.getText());
     }
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

  public void setAnaObjStack(Stack<WlzObject> stack){
    _anaObjStack = stack;
  }

  public void setAnaFilesStack(Stack<File> stack){
    _anaFileStack = stack;
  }

  public void setPathLengthToAnatomy(int length){
    _pathLengthToAnatomy = length;
  }

//----------------------------------------------------
  public Stack<File> getAnaFileStack() { return _anaFileStack; }
  public Stack<File> getEmbFileStack() { return _embFileStack; }
  public Stack<File> getXEmbFileStack() { return _xembFileStack; }
  public Stack<WlzObject> getAnaObjStack() { return _anaObjStack; }
  public Stack<WlzObject> getEmbObjStack() { return _embObjStack; }
  public Stack<WlzObject> getXEmbObjStack() { return _xembObjStack; }

  public SelectableMenu getAnaMenu() { return _anaMenu; }
  public SelectableMenu getEmbMenu() { return _embMenu; }
  public SelectableMenu getXEmbMenu() { return _xembMenu; }

  public int getPathLengthToAnatomy() { return _pathLengthToAnatomy; }

} // class AnatomyBuilder
