/* file name  : NodeImageContainer.java
 * created    : Mon 03 Sep 2007 09:35:42 BST
 */
package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;
import java.util.*;
import javax.imageio.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/** 
 * Image container subclass that displays node (heatmap) images.
 * 
 */
public class NodeImageContainer extends ImageContainer {
   String movieURI;

   public NodeImageContainer(String s, BufferedImage i, String c, int sortIndex, String movieURI) {

      super(s,i,sortIndex);
      correlation = c;
      this.movieURI = movieURI;

      //System.out.println("HeatmapImageContainer: BufferedImage is "+getWidth()+" x "+getHeight());

      labels.add(new ClickableLabel(nodeID));
      if (ImageViewManager.isDebugMode()) {
	 labels.add(new ClickableLabel("Sort index: " + sortIndex));
      }
      if (correlation != null) {
	 labels.add(new ClickableLabel("Correlation: " + correlation));
      }
      addLabels();
      this.setBackground(Color.orange);
   }

   public void paintComponent(Graphics g) {
      super.paintComponent(g);

      if (imagePanel != null) {
	 System.out.println("HeatmapImageContainer: paintComponent selected = "+selected);
	 if (selected) {
	    imagePanel.setBorder(BorderFactory.createLineBorder(Color.red,borderWidth));
	 } else {
	    imagePanel.setBorder(BorderFactory.createLineBorder(Color.green,1));
	    //imagePanel.setBorder(BorderFactory.createEmptyBorder());
	 }
      }
      //System.out.println("HeatmapImageContainer: exit paintComponent()");
   }

   public String getHyperlinkTarget() {
      return movieURI;
   }
}
