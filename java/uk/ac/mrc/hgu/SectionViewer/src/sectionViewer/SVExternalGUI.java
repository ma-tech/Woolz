package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

// import hguUntil.*;

public class SVExternalGUI extends JFrame {

  private final boolean _debug = false;

  JPanel _contentPane;
  JRootPane _rootPane;

  // these need to be visible outside of guiInit()
  int totalH = 500;
  int totalW = 400;
  int pad = 1;
  //=========================================================
  // constructor
  //=========================================================
  public SVExternalGUI(String str) {

    super(str);

    if(_debug == true) System.out.println("enter SVExternalGUI");

    try {
      guiInit();
    }
    catch(Exception e) {
      System.out.println("SVExternalGUI");
      e.printStackTrace();
    }

    if(_debug) System.out.println("exit SVExternalGUI");
  }

//========================================
  private void guiInit() throws Exception {

    if(_debug) System.out.println("enter guiInit");

    _contentPane = (JPanel)this.getContentPane();
    _rootPane = this.getRootPane();

    _contentPane.setLayout(new BorderLayout());
    _contentPane.setPreferredSize(new Dimension(totalW+pad, totalH+pad));
    _contentPane.setMinimumSize(new Dimension(totalW+pad, totalH+pad));

    //----------------------------------------------------------
    if(_debug) System.out.println("exit guiInit");
  } // guiInit()


//----------------------------------------------------------

} // class SVExternalGUI
