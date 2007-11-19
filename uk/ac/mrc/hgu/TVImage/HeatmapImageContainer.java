/* file name  : HeatmapImageContainer.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
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
 * Image container subclass that displays heatmap images.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
public class HeatmapImageContainer extends ImageContainer {

   String movieURI;

//-------------------------------------------------------------------
   public HeatmapImageContainer(String s, BufferedImage i, String c, int sortIndex, String movieURI) {

      super(s,i,sortIndex);
      correlation = c;
      this.movieURI = movieURI;

      labels.add(new ClickableLabel(nodeID));
      if (ImageViewManager.isDebugMode()) {
	 labels.add(new ClickableLabel("Sort index: " + sortIndex));
      }
      if (correlation != null) {
	 labels.add(new ClickableLabel("Correlation: " + correlation));
      }
      addLabels();
   }
//-------------------------------------------------------------------
   public String getHyperlinkTarget() {
      if(movieURI == null || movieURI.equals("")) {
         movieURI = "http:www.hgu.mrc.ac.uk";
      }
      return movieURI;
   }
//-------------------------------------------------------------------
}
