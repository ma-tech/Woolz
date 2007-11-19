package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.util.*;
import javax.imageio.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/** 
 * JScrollPane that holds a single ImageContainer.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
class TopImageScrollPane extends ImageScrollPane {

	String nodename = "";

	/** 
	 * Creates a new JScrollPane to display images in a
	 * GridLayout, and sets preferred sizes appropriately.
	 * 
	 * @param group the group of synchronised scrolling areas that
	 * this scrollpane should become a member of
	 * @param type the type of images that this scroll pane should
	 * display
	 */
	TopImageScrollPane(ScrollGroup group, int type, ImageServer imgProducer, ImageViewManager ivm) {

		super(group, type, imgProducer, ivm);
		//content.setBackground(Color.green);

	} // TopImageScrollPane


	void setPanelSizes() {
	   if (panels.size() < 1) {
	      return;
	   }

	   // SPECIAL CASE: for single node selections
	   if (panels.size() == 1) {
	      if (type == ImageViewManager.TYPE_NODE ||
		    type == ImageViewManager.TYPE_RAW ||
		    type == ImageViewManager.TYPE_ANNOTATED) {
		 for (ImageContainer c : panels) {
		    c.fitToSize(new Dimension(300, 300));
		 }
		 panelSize = panels.iterator().next().getPreferredSize();
	      }
	   } else {
	      for (ImageContainer c : panels) {
	         if(type == ImageViewManager.TYPE_NODE) {
		    c.fitToSize(new Dimension(200,200));
		 } else {
		    c.fitToSize(ImageView.MAX_PANEL_SIZE);
		 }
	      }
	      panelSize = ImageView.MAX_PANEL_SIZE;
	   }
	} // setPanelSizes

} // TopImageScrollPane
