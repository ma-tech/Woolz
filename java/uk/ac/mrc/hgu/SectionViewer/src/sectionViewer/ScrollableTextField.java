package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

/**
 *   A text field which can be scrolled even if it is non-editable.
 *   Scrolling is achieved by dragging the mouse left or right
 *   inside the text field.
 */
public class ScrollableTextField extends JTextField {

  /**   true if text is scrolled to the left */
  boolean scrollLeft = false;

  /**   the position of the mouse inside the text field. */
  int X;

  /**   the previous position of the mouse inside the text field. */
  int XPrev = 0;

  /**   the number of pixels to scroll the text each time. */
  final int  DELTA = 5;

  /**
   *   Creates a ScrollableTextField with no initial text.
   */
  public ScrollableTextField() {
    super();
  }

  /**
   *   Creates a ScrollableTextField with the given initial text.
   *   @param str the initial text.
   */
  public ScrollableTextField(String str) {
    super(str);
  }

  /**
   *   Creates a ScrollableTextField with the given number of columns.
   *   @param clms the initial number of columns.
   */
  public ScrollableTextField(int clms) {
    super(clms);
  }

  /**
   *   Creates a ScrollableTextField with
   *   the given initial text and
   *   the given number of columns.
   *   @param str the initial text.
   *   @param clms the initial number of columns.
   */
  public ScrollableTextField(String str, int clms) {
    super(str, clms);
  }

  /**
   *   Event Handler for mouse movement in the text field.
   */
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

  /**
   *   Sets the scrollOffset for the text field.
   *   Text will scroll left or right by the value of DELTA.
   */
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

