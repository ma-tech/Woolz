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

public class SVInternal
    extends SVInternalGUI {

  private final boolean _debug = false;

  protected SVPanel _svPanel = null;


  //=========================================================
  // constructor
  //=========================================================
  public SVInternal(String viewstr) {

    super(viewstr);

    if (_debug) System.out.println("enter SVInternal");

    if (_debug)
      System.out.println("exit SVInternal");
  } // constructor

  public SVPanel getSVPanel() {
     return _svPanel;
  }

  public void setSVPanel(SVPanel svp) {
     _svPanel = svp;
     _contentPane.add(_svPanel, BorderLayout.CENTER);
  }


} // class SVInternal
