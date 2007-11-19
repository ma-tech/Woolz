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
 * JScrollPane that holds a collection of ImageContainers.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
class ImageScrollPane extends JScrollPane implements Observer {

   JPanel content = null;
   SortedSet<ImageContainer> panels = new TreeSet<ImageContainer>(new PanelComparator());
   Dimension panelSize;
   GridLayout splm = null;
   Dimension available;
   ImageServer imgProducer;
   ImageView parent;
   int type = ImageViewManager.TYPE_NULL;
   String nodename = "";
   private boolean _debug = false;
   //-----------------------------------------------------------------------------------------
   /** 
    * Creates a new JScrollPane to display images in a
    * GridLayout, and sets preferred sizes appropriately.
    * 
    * @param group the group of synchronised scrolling areas that
    * this scrollpane should become a member of
    * @param type the type of images that this scroll pane should
    * display
    */
   ImageScrollPane(ScrollGroup group, int type, ImageServer imgProducer, ImageViewManager ivm) {

      super();
      content = new JPanel();
      //content.setBackground(Color.red);
      //content.setBorder(BorderFactory.createLineBorder(Color.blue,1));
      content.setBackground(ImageViewManager.BGCOLOR);
      setViewportView(content);
      addMouseWheelListener(group);
      this.type = type;
      this.imgProducer = imgProducer;
      if (group != null) {
	 JScrollBar jsb = getVerticalScrollBar();
	 group.put(new Integer(type), jsb);
	 jsb.addAdjustmentListener(group);
	 jsb.addMouseListener(group);
	 jsb.addMouseWheelListener(group);
      }
   } // ImageScrollPane
   //-----------------------------------------------------------------------------------------
   public void setParent(ImageView iv) { parent = iv; }
   //-----------------------------------------------------------------------------------------
   public boolean isEmpty() { return panels.isEmpty(); }
   //-----------------------------------------------------------------------------------------
   protected String getType() {
      if(type == ImageViewManager.TYPE_NODE) {
	 return "Node";
      } else if(type == ImageViewManager.TYPE_HEATMAP) {
	 return "Heatmap";
      } else if(type == ImageViewManager.TYPE_RAW) {
	 return "Raw";
      } else if(type == ImageViewManager.TYPE_ANNOTATED) {
	 return "Annotated";
      } else if(type == ImageViewManager.TYPE_NULL) {
	 return "Null";
      }
      return "not known";
   }
   //-----------------------------------------------------------------------------------------
   public SortedSet<ImageContainer> getPanels() {
      return panels;
   }
   //-----------------------------------------------------------------------------------------
   /** 
    * Returns a list of node identifiers for the displayed images
    */
   public SortedSet<String> getIDs() {
      SortedSet<String> ids = new TreeSet<String>();
      for (ImageContainer c : panels) {
	 ids.add(c.getID());
      }
      return ids;
   }
   //-----------------------------------------------------------------------------------------
   /** 
    * Prepares a new ImageContainer and adds it to a list of
    * all the ImageContainers that are waiting to be drawn.
    *
    * @param id The name of the node to get (GENE??X or NODE??X)
    */
   void addPanel(String id) {
      assert id != null;

      if (id == null) {
	 (new IllegalArgumentException("Tried to add a null image panel!")).printStackTrace();
      }
      //System.out.println("ImageScrollPane: enter addPanel for "+getType()+", "+id);

      ImageContainer container = null;
      container = imgProducer.getMedia(id);
      assert container != null;

      container.setParent(parent);
      panels.add(container);

      setStuff();
      //System.out.println("ImageScrollPane: exit addPanel for "+getType()+", "+id);
   }
   //-----------------------------------------------------------------------------------------
   void removePanels(Collection<String> ids) {
      assert ids != null;
      for (String s : ids) {
	 removePanel(s);
      }

      setStuff();
   }
   //-----------------------------------------------------------------------------------------
   void removePanel(String id) {
      assert id != null;
      if (id == null) {
	 return;
      }

      ImageContainer temp = null;

      for (ImageContainer p : panels) {
	 if (p.getID().equals(id)) {
	    temp = p;
	 }
      }
      if (temp != null) {
	 panels.remove(temp);
	 content.remove(temp);
      }
   }
   //-----------------------------------------------------------------------------------------
   protected void clear() {
      content.removeAll();
      panels.clear();
   }
   //-----------------------------------------------------------------------------------------
   void setStuff() {
      doLocalLayout(null);
      setViewportView(content);
      revalidate();
      repaint();
   }

   //-----------------------------------------------------------------------------------------
   /**
    *  Calculates number of columns for the grid.
    *  Don't set the number of rows or columns will be ignored.
    */
   void doLocalLayout(Dimension available) {

      //System.out.println(">>>>>> ImageScrollPane: enter doLocalLayout");
      if(available != null) {
	 this.available = available;
	 //System.out.println("available: "+available.width+", "+available.height);
      }
      int cols = 4;
      int hgap = 2;
      int vgap = 2;
      int totalImgs = panels.size();

      if (panels.isEmpty()) {
	 setPreferredSize(available);
	 content.setPreferredSize(available);
      } else {
	 cols = calcCols();
	 splm = new GridLayout(0, cols);
	 splm.setHgap(hgap);
	 splm.setVgap(vgap);
	 content.setLayout(splm);
	 content.removeAll();
	 //System.out.println("ImgSP: "+getType()+" doLocalLayout: content.removeAll() ---------------------------------------------------");
	 for (ImageContainer c : panels) {
	    //System.out.println("ImgSP: "+getType()+" doLocalLayout: content.add() "+c.getID()+" -------------------------------------------");
	    content.add(c);
	 }
      }
      //System.out.println("<<<<<< ImageScrollPane: exit doLocalLayout");
   }
   //-----------------------------------------------------------------------------------------
   private int calcCols() {

      int cols = 4;
      ImageContainer ic = null;

      if(panels == null || panels.isEmpty()) {
	 return cols;
      }
      ic = panels.first();
      if(ic == null) {
	 return cols;
      }
      int icWidth = ic.getImageContainerWidth();
      //System.out.println("Image container width = "+icWidth);

      if(available == null) {
	 return cols;
      }
      int thisWidth = available.width;
      //System.out.println("Image scrollPane width = "+thisWidth);

      cols = thisWidth / icWidth;
      //System.out.println("cols = "+cols);

      return cols;
   }
   //-----------------------------------------------------------------------------------------
   public void update(Observable o, Object arg) {
      //System.out.println(">>>>>> ImageScrollPane: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
      //System.out.println("<<<<<< ImageScrollPane: exit update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
   }
   //-----------------------------------------------------------------------------------------
   public JPanel getContent() {
      return content;
   }



} // ImageScrollPane
