package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

public class SVInternalGUI extends JInternalFrame {

  private final boolean _debug = false;

  JPanel _contentPane;

  // these need to be visible outside of guiInit()
  int totalH = 500;
  int totalW = 400;
  int pad = 1;
  //=========================================================
  // constructor
  //=========================================================
  public SVInternalGUI(String str) {

    super(str);

    if(_debug == true) System.out.println("enter SVInternalGUI");

    try {
      guiInit();
    }
    catch(Exception e) {
      System.out.println("SVInternalGUI");
      e.printStackTrace();
    }

    if(_debug) System.out.println("exit SVInternalGUI");
  }

//========================================
  private void guiInit() throws Exception {

    if(_debug) System.out.println("enter guiInit");

    _contentPane = (JPanel)this.getContentPane();

    _contentPane.setLayout(new BorderLayout());
    _contentPane.setPreferredSize(new Dimension(totalW+pad, totalH+pad));
    _contentPane.setMinimumSize(new Dimension(totalW+pad, totalH+pad));

    //----------------------------------------------------------
    if(_debug) System.out.println("exit guiInit");
  } // guiInit()


//----------------------------------------------------------

} // class SVInternalGUI
