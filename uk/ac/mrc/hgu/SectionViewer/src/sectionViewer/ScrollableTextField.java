package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

public class ScrollableTextField extends JTextField {

  boolean scrollLeft = false;
  int X;
  int XPrev = 0;
  final int  DELTA = 5;

   /* constructors */
  public ScrollableTextField() {
    super();
  }
  public ScrollableTextField(String str) {
    super(str);
  }
  public ScrollableTextField(int clms) {
    super(clms);
  }
  public ScrollableTextField(String str, int clms) {
    super(str, clms);
  }

  public void processMouseMotionEvent (MouseEvent e) {
    switch(e.getID()) {
      case MouseEvent.MOUSE_DRAGGED:
        X = e.getX();
        if(X > XPrev) {
          scrollLeft = false;
        } else {
          scrollLeft = true;
        }
        XPrev = X;
        scroll();
        break;
      default:
        break;
    }
  }

  public void scroll() {
    int offset = super.getScrollOffset();
    //System.out.println("scrollOffset = "+ offset);
    if(scrollLeft == false) {
      if(offset < DELTA) {
        offset = 0;
      } else {
        offset = offset - DELTA;
      }
    } else {
      offset = offset + DELTA;
    }
    //System.out.println("setting scrollOffset = "+ offset);
    super.setScrollOffset(offset);
  }

}

