package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
/*
import java.awt.image.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.imageio.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;
*/

/** 
 *
 * A panel that directly contains images to display.
 * 
 */
public class ImagePanel extends JPanel implements MouseListener, TVTypes {

   private int borderWidth = 2;
   private ImageContainer parent;
   private Color bgc = TVTypes.BGCOLOR;
   protected java.util.Timer tim = null;
   protected boolean doubleclick = false;
   private int doubleClickTime = 500; // milliseconds

   BufferedImage img = null;
   Dimension actualImageSize = new Dimension(100, 100);

   public ImagePanel() {
      addMouseListener(this);
   }

   public void paintComponent(Graphics g) {
      assert img != null;
      assert parent != null;
      Dimension siz = parent.getImagePanelSize();
      Color col = g.getColor();
      g.setColor(bgc);
      g.fillRect(0,0,siz.width,siz.height);
      g.setColor(col);

      Point pnt = getOrigin();
      if(img != null) {
	 g.drawImage(img, pnt.x, pnt.y, actualImageSize.width, actualImageSize.height, bgc, this);
      }
      //System.out.println("ImagePanel: ID = "+parent.getID()+" exit paintComponent()");
   }

   public void setImage(BufferedImage bi) {
      if(parent == null) {
         //System.out.println("ImagePanel.setImage: returning, parent is null");
         return;
      }
      if(bi == null) {
         //System.out.println("ImagePanel.setImage: returning, img is null");
         return;
      }
      img = bi;
      actualImageSize = new Dimension(img.getWidth(), img.getHeight());
   }

   public BufferedImage getImage() {
      return img;
   }

   public void setParent(ImageContainer p) {
      parent = p;
      bgc = p.getBGC();
   }

   public ImageContainer getParent() {
      return parent;
   }

   public void setBorderWidth(int w) {
      borderWidth = w;
   }

   public Dimension getMaximumSize() {
      //System.out.println("getMaximumSize: "+getPreferredSize().toString());
      return getPreferredSize();
   }

   public Dimension getMinimumSize() {
      //System.out.println("getMinimumSize: "+getPreferredSize().toString());
      return getPreferredSize();
   }

   /**
    *   We need to differentiate between single and double mouse clicks.
    *   If the user double-clicks we don't want to do the action for single-click first.
    *   A timer schedules single click unless a second click arrives within 'doubleClickTime' mS.
    */
   public void mouseClicked(MouseEvent me) {
      tim = new java.util.Timer();
      if (me.getClickCount() >= 2){
         doubleclick = true;
      } else {
         doubleclick = false;
	 tim.schedule(new MouseClickTask(), doubleClickTime);
      }
   }
   public void mousePressed(MouseEvent me) {
      //System.out.println("mouse pressed");
   }
   public void mouseReleased(MouseEvent me) {
      //System.out.println("mouse released");
   }
   public void mouseEntered(MouseEvent me) {
      String path = parent.getHyperlinkTarget();
      if (path != null && path != "") {
	 Cursor cursor = getCursor();
	 setCursor(cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
      }
   }

   public void mouseExited(MouseEvent me) {
      String path = parent.getHyperlinkTarget();
      if (path != null && path != "") {
	 Cursor cursor = getCursor();
	 setCursor(cursor.getDefaultCursor());
      }
   }

   private Point getOrigin() {

      Dimension dim = null;

      if(parent == null) {
         return new Point(borderWidth, borderWidth);
      } else {
	 dim = parent.getImagePanelSize();
      }
      double xcent = dim.getWidth() / 2.0;
      double ycent = dim.getHeight() / 2.0;

      int xorig = (int)(xcent - (actualImageSize.getWidth() / 2.0));
      int yorig = (int)(ycent - (actualImageSize.getHeight() / 2.0));

      return new Point(xorig, yorig);
   }
//---------------------------------------------------------------------------
   class MouseClickTask extends TimerTask{
      public void run() {
         if(doubleclick) {
	    parent.doDoubleMouseClick();
	    //System.out.println("double click");
	 } else {
	    parent.doSingleMouseClick();
	    //System.out.println("single click");
	 }
	 tim.cancel();
      }
   }

} // class ImagePanel
