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

/** 
 * Superclass for image panels displayed in an ImageScrollPane.
 * Consists of a JPanel that holds another JPanel for the actual
 * image, a collection of labels and a collection of other stuff
 * (generally empty).
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
public abstract class ImageContainer extends JPanel {

	protected int borderWidth = 2;
	protected BufferedImage img = null;
	protected ImagePanel imagePanel = null;
	protected JPanel labelPanel = null;
	protected ArrayList<ClickableLabel> labels = new ArrayList<ClickableLabel>();
	protected ArrayList<JComponent> otherStuff = new ArrayList<JComponent>();
	protected String nodeID = "";
	protected String correlation = null;
	protected boolean selected = false;
	protected ImageView tviv;
	private int sortIndex = 0;
	private int thisW = 140;
	private int thisH = 185;
	private int imagePanelW = 140;
	private int imagePanelH = 140;
	private int labelPanelW = 140;
	private int labelPanelH = 35;
	private Color selectedCol = Color.red;
	private Color notSelectedCol = Color.green;
	//private Color bgc = new Color(240,240,245);
	private Color bgc = ImageViewManager.BGCOLOR;
	private boolean _debug = false;

//---------------------------------------------------------------------------
	public int getSortIndex() { return sortIndex; }
//---------------------------------------------------------------------------
	/**
	 * In the constructor of your subclass, call this first to set
	 * some defaults.  Then call finishLoading() to make sure your
	 * own settings are applied.
	 *
	 * @param nodeID node ID
	 * @param i image to be displayed
	 * @param maxWidth
	 */
	public ImageContainer(String nodeID, BufferedImage i, int sortIndex) {
	   this.nodeID = nodeID;
	   this.sortIndex = sortIndex;
	   this.img = i;
	   if (i == null) {
	      System.out.println("ImageContainer "+nodeID+" could not load image");
	      labels.add(new ClickableLabel("Could not load image"));
	   }
	   initPanels();
	}
//---------------------------------------------------------------------------
	private void initPanels() {
	   this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
	   this.setBackground(ImageViewManager.BGCOLOR);
	   //this.setBorder(BorderFactory.createLineBorder(Color.blue, borderWidth));
	   this.setPreferredSize(new Dimension(thisW, thisH));
	   this.setMinimumSize(new Dimension(thisW, thisH));
	   this.setMaximumSize(new Dimension(thisW, thisH));

	   imagePanel = new ImagePanel();
	   imagePanel.setParent(this);
	   imagePanel.setPreferredSize(new Dimension(imagePanelW, imagePanelH));
	   imagePanel.setMaximumSize(new Dimension(imagePanelW, imagePanelH));
	   imagePanel.setMinimumSize(new Dimension(imagePanelW, imagePanelH));
	   imagePanel.setImage(img);
	   imagePanel.setBorderWidth(borderWidth);
	   imagePanel.setAlignmentX(Component.CENTER_ALIGNMENT);
	   //imagePanel.setBorder(BorderFactory.createLineBorder(notSelectedCol, borderWidth));
	   //imagePanel.setBackground(ImageViewManager.BGCOLOR);
	   imagePanel.setBackground(bgc);

	   labelPanel = new JPanel();
	   labelPanel.setLayout(new BoxLayout(labelPanel, BoxLayout.Y_AXIS));
	   labelPanel.setPreferredSize(new Dimension(labelPanelW, labelPanelH));
	   labelPanel.setMaximumSize(new Dimension(labelPanelW, labelPanelH));
	   labelPanel.setMinimumSize(new Dimension(labelPanelW, labelPanelH));
	   //labelPanel.setBackground(ImageViewManager.BGCOLOR);
	   labelPanel.setBackground(bgc);

	   add(Box.createRigidArea(new Dimension(0,2)));
	   add(imagePanel);
	   add(Box.createRigidArea(new Dimension(0,2)));
	   add(labelPanel);
	}

//---------------------------------------------------------------------------
	protected void addLabels() {
	   for (ClickableLabel il : labels) {
	      il.setAlignmentX(Component.CENTER_ALIGNMENT);
	      labelPanel.add(il);
	   }
	}
//---------------------------------------------------------------------------
	public boolean equals(ImageContainer c) {
	   return (this.getID().equals(c.getID()));
	}
//---------------------------------------------------------------------------
	/** 
	 * Resizes the image container to fit within the specified
	 * dimension.
	 * 
	 * @param d dimension that the image container must fit within
	 */
	public void fitToSize(Dimension d) {
	   if (d == null || d.width < 0 || d.height < 0) {
	      (new IllegalArgumentException()).printStackTrace();
	   }

	   assert imagePanel != null;

	   int labelHeight = 0;
	   for (ClickableLabel il : labels) {
	      il.fitToSize(d);
	      labelHeight += il.getPreferredSize().height;
	   }
	   for (JComponent jc : otherStuff) {
	      labelHeight += jc.getPreferredSize().height;
	   }

	   if (imagePanel != null) {
	      //imagePanel.fitToSize( new Dimension(d.width, d.height-labelHeight));
	   }
	   revalidate();
	   repaint();
	}
//---------------------------------------------------------------------------
	public void setParent(ImageView v) {
	   tviv = v;
	   //System.out.println("ImageContainer: "+getID()+" setParent");
	}
//---------------------------------------------------------------------------
	public String getID() { return nodeID; }
//---------------------------------------------------------------------------
	public String toString() { return getID(); }
//---------------------------------------------------------------------------
	public BufferedImage getImage() { return imagePanel.getImage(); }
//---------------------------------------------------------------------------
	public boolean isSelected() { return selected; }
//---------------------------------------------------------------------------
	public void setSelected(boolean b) {
	   if(_debug) {
	      System.out.println("ImageContainer: ID = "+getID()+" setSelected "+b);
	   }
	   selected = b;
	   if(selected) {
	      imagePanel.setBorder(BorderFactory.createLineBorder(selectedCol, borderWidth));
	   } else {
	      imagePanel.setBorder(null);
	   }
	   imagePanel.invalidate();
	   imagePanel.revalidate();
	}
//---------------------------------------------------------------------------
	abstract String getHyperlinkTarget();
//---------------------------------------------------------------------------
	protected int getImageContainerWidth() {
	   return thisW;
      }
//---------------------------------------------------------------------------
   /**
    *   Makes a scaled BufferedImage.
    *   (based on comp.lang.java.gui reply by Shannon Hickey - Swing Team)
    *   Google:	scaling a BufferedImage w/ Java 2D
    *   @param bi the unscaled image.
    *   @return the scaled image.
    */
   protected BufferedImage scaleImg(BufferedImage bi) {

      int H = bi.getHeight();
      int W = bi.getWidth();
      //System.out.println("ImageContainer.scaleImg: panel is "+imagePanelW+" wide X "+imagePanelH+" high");
      //System.out.println("ImageContainer.scaleImg: img is "+W+" wide X "+H+" high");

      if(H <= imagePanelH && W <= imagePanelW) {
	 // no need to scale image to fit container
	 //System.out.println("ImageContainer.scaleImg: returning, img fits container already");
         return bi;
      }

      double ratio = ((double)H/(double)W);
      boolean portrait = (ratio >= 1.0); // the image is taller than it is wide (or square).
      //System.out.println("Image H/W = "+ratio);

      int scaledW = 0;
      int scaledH = 0;
      double newRatio = 1.0;

      if(portrait) {
	 //System.out.println("portrait");
	 scaledH = imagePanelH;
	 scaledW = (int)(Math.ceil((imagePanelH / ratio)));
	 //System.out.println("scaled img is "+scaledW+" wide X "+scaledH+" high");
	 newRatio = ((double)scaledW/(double)imagePanelW);
	 //System.out.println("new ratio = "+ratio);
	 if(newRatio > 1.0) {
	    scaledH /= newRatio;
	    scaledW /= newRatio;
	 }
	 //System.out.println("scaled img is now "+scaledW+" wide X "+scaledH+" high");
      } else {
	 //System.out.println("landscape");
	 scaledW = imagePanelW;
	 scaledH = (int)(Math.ceil((imagePanelW * ratio)));
	 //System.out.println("scaled img is "+scaledW+" wide X "+scaledH+" high");
	 newRatio = ((double)scaledH/(double)imagePanelH);
	 //System.out.println("new ratio = "+ratio);
	 if(newRatio > 1.0) {
	    scaledH /= newRatio;
	    scaledW /= newRatio;
	 }
	 //System.out.println("scaled img is now "+scaledW+" wide X "+scaledH+" high");
      }

      BufferedImage scaledImg = new BufferedImage(scaledW, scaledH, bi.getType());
      Graphics2D g = (Graphics2D)(scaledImg.getGraphics());

      // Frame implements imageObserver interface.
      Frame dummyObserver = new Frame();

      //g.setBackground(Color.magenta);
      g.clearRect(0, 0, scaledW, scaledH);
      g.drawImage(bi, 0, 0, scaledW, scaledH, dummyObserver);
      g.dispose();

      return scaledImg;
   } // scaleImg()
//---------------------------------------------------------------------------
   /**
    *   Handles single mouse click events from images.
    *   If the correlation slider (or text box) has been changed:
    *      single-clicking a heatmap image (Node) highlights the corresponding part of the tree.
    *      single-clicking a heatmap image (Gene) does nothing to tree.
    *   If the tree has been clicked:
    *      ???
    */
   protected void doSingleMouseClick() {
      if(_debug) {
	 System.out.println(">>>>>> ImageContainer: "+getID()+" enter doSingleMouseClick");
      }
      boolean wasSelected = selected;
      if(tviv.getManager().isCorrelationView()) {
	 //System.out.println("isCorrelationValueSet() "+tviv.getManager().getGeneSelection().isCorrelationValueSet());
	 //System.out.println("getSelectedNodes().size() "+tviv.getManager().getGeneSelection().getSelectedNodes().size());
	 int index = -1;
	 tviv.setSelectedContainer(null);
	 if(!wasSelected) {
	    tviv.setSelectedContainer(this);
	 }
      }
      if(_debug) {
	 System.out.println("<<<<<< ImageContainer: "+getID()+" exit doSingleMouseClick");
      }
   }
//---------------------------------------------------------------------------
   /**
    *   Handles double (or more) mouse click events from images.
    *   If the correlation slider (or text box) has been changed:
    *      double-clicking a heatmap image (Node) corresponds to clicking that node on the tree.
    *   If the tree has been clicked:
    *      ???
    */
   protected void doDoubleMouseClick() {
      if(_debug) {
	 System.out.println(">>>>>> ImageContainer: "+getID()+" enter doDoubleMouseClick");
      }
      String path = "";
      if(isGene(getID())) {
	 tviv.getManager().setDoubleClickGene(true);
      } else {
	 tviv.getManager().setDoubleClickNode(true);
      }
      //if (tviv.getManager().getGeneSelection().isCorrelationValueSet() ||
         //tviv.getManager().getGeneSelection().getSelectedNodes().size() > 1) {
      if(tviv.getManager().isCorrelationView()) {
	 //System.out.println("isCorrelationValueSet() "+tviv.getManager().getGeneSelection().isCorrelationValueSet());
	 //System.out.println("getSelectedNodes().size() "+tviv.getManager().getGeneSelection().getSelectedNodes().size());
	 tviv.getManager().getGeneSelection().resetCorrelationValue();
	 int index = -1;
	 tviv.setSelectedContainer(null);
	 //tviv.getManager().getGeneSelection().setSelectedNode(getID());
	 if (isGene(getID())) {
	    tviv.getManager().getGeneSelection().setSelectedNode(null);
	    index = tviv.getManager().getGeneIndex(getID());
	    tviv.getManager().getGeneSelection().selectIndexRange(index,index);
	 } else {
	    tviv.getManager().getGeneSelection().setSelectedNode(getID());
	 }
	 tviv.getManager().getGeneSelection().notifyObservers();
      } else {
	 path = getHyperlinkTarget();
	 if (path != null && path != "") {
	    System.out.println(path);
	    tviv.getManager().displayURL(path);
	 }
      }
      tviv.getManager().setDoubleClick(false);
      if(_debug) {
	 System.out.println("<<<<<< ImageContainer: "+getID()+" exit doDoubleMouseClick");
      }
   }
//---------------------------------------------------------------------------
   protected boolean isGene(String id) {
      boolean ret = false;
      String str = id.toLowerCase();
      if(str.startsWith("gene")) {
         ret = true;
      }
      return ret;
   }
//---------------------------------------------------------------------------
   protected Dimension getImagePanelSize() {
      return new Dimension(imagePanelW, imagePanelH);
   }
//---------------------------------------------------------------------------
   protected Color getBGC() {
      return bgc;
   }
//---------------------------------------------------------------------------

}
