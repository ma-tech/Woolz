/* file name  : ClickableLabel.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
 * created    : Mon 03 Sep 2007 09:13:18 BST
 */
package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.imageio.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

import external.*;

/** 
 * Extends the JLabel class to provide hyperlink-style
 * functionality.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
class ClickableLabel extends JLabel implements MouseListener
{
   Color defaultColor = Color.black;
   Color highlightedColor = Color.blue;
   String url;
   ImageViewManager ivm;

   /** 
    * Create a new ClickableLabel.  If the specified URL is null,
    * then the label should not be clickable.
    * 
    * @param s string to be displayed
    * @param url URL to be opened when clicked
    * @param ivm parent image view manager, required for its link
    * to the displayURL method in
    * edu.stanford.genetics.treeview.{Tree,}ViewFrame.java
    */
   ClickableLabel(String s, String url, ImageViewManager ivm) {
      super(s);
      this.url = url;
      this.ivm = ivm;
      if (url != null) {
	 addMouseListener(this);
	 //system.out.println("ClickableLabel: "+s+", url = "+url);
      }
   }

   /** 
    * Create a non-clickable ClickableLabel
    * 
    * @param s string to display
    */
   ClickableLabel(String s) { 
      this(s,null,null);
   }

   public void paintComponent(Graphics g) {
      Graphics2D g2 = (Graphics2D)g;
      g2.setRenderingHint(
	    RenderingHints.KEY_TEXT_ANTIALIASING,
	    RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
      super.paintComponent(g2);
   }

   public void fitToSize(Dimension d) {
      while (getPreferredSize().width > d.width) {
	 Font oldFont = getFont();
	 setFont(new Font(
		  oldFont.getName(),
		  oldFont.getStyle(),
		  oldFont.getSize() - 1));
      }
   }

   public void mouseClicked(MouseEvent e) {
      ivm.displayURL(url);
      /*
      try {
	 BrowserLauncher.openURL(url);
	 //BrowserControl.getBrowserControl().displayURL(url);
      }
      catch (MalformedURLException mue) {
	 mue.printStackTrace();
      }
      catch (IOException ioe) {
	 ioe.printStackTrace();
      }
      */
   }

   public void mousePressed(MouseEvent e) {}
   public void mouseReleased(MouseEvent e) {}

   public void mouseEntered(MouseEvent e) {
      setForeground(highlightedColor);
      Cursor cursor = getCursor();
      setCursor(cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
   }

   public void mouseExited(MouseEvent e) {
      setForeground(defaultColor);
      Cursor cursor = getCursor();
      setCursor(cursor.getDefaultCursor());
   }
}
