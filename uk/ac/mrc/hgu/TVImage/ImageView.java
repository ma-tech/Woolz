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
 * Displays images appropriate to the selected gene(s), such as
 * assay images and/or heatmaps.  Generally the application will
 * instantiate a few of these, each of which loads a different
 * type of image (heatmap/raw/annotated...)
 *
 * Since this class is a JSplitPane, it may hold two
 * subcomponents.  These are currently assumed to be
 * ImageScrollPanes.
 *
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 *
 * Changes
 * Paul Smith - Change to allow for raw and annotated images
 */
class ImageView extends JSplitPane implements Observer, ComponentListener, TVTypes {

   private ImageViewManager tvim;
   private Set<String> topPanels = new TreeSet<String>();
   private Set<String> bottomPanels = new TreeSet<String>();
   private ImageContainer selectedImageContainer = null;
   public static final Dimension MAX_PANEL_SIZE = new Dimension(120,180);
   private String description = "";
   private ImageServer imgProducer = null;
   public static final int FUDGE_FACTOR = 10;
   public static final boolean INCREMENTAL_IMAGE_LOADING = true;
   private int defaultDividerSize;
   private ImageScrollPane top;
   private ImageScrollPane bottom;
   private boolean displaysSelectionByCorrelation;
   private boolean _debug = false;

//-------------------------------------------------------------------
   public ImageView(
	 ImageViewManager m,
	 XmlAppSettings appSettings,
	 boolean c,
	 String d,
	 ImageScrollPane top,
	 ImageScrollPane bottom) {

      super(JSplitPane.VERTICAL_SPLIT);

      assert m != null;
      assert d != null;

      //System.out.println("ImageView: constructor, c = "+c+", description = "+d);

      this.tvim = m;
      this.displaysSelectionByCorrelation = c;
      this.description = d;
      this.top = top;
      this.bottom = bottom;

      defaultDividerSize = getDividerSize();
      setOneTouchExpandable(true);

      addComponentListener(this);

      top.setParent(this);
      bottom.setParent(this);
      setTopComponent(top);
      setBottomComponent(bottom);
   }

//-------------------------------------------------------------------
   public void componentResized(ComponentEvent e) {
      //System.out.println("ImageView now "+getSize().toString());
      top.doLocalLayout(getSize());
      bottom.doLocalLayout(getSize());
      doMediaLayout();
   }
   //-------------------------------------------------------------------
   public void componentHidden(ComponentEvent e) {}
   //-------------------------------------------------------------------
   public void componentMoved(ComponentEvent e) {}
   //-------------------------------------------------------------------
   public void componentShown(ComponentEvent e) {}
   //-------------------------------------------------------------------
   public String viewName() {
      return "ImageView";
   }
   //-------------------------------------------------------------------
   public String viewDescription() {
      return description;
   }
   //-------------------------------------------------------------------
   public ImageViewManager getManager() {
      return tvim; 
   }
   //-------------------------------------------------------------------
   /** 
    * Returns the image container that is currently selected.
    * 
    * @return the currently selected image container
    */
   public ImageContainer getSelectedContainer() {
      return selectedImageContainer;
   }
   //-------------------------------------------------------------------
   /** 
    * Selects a ImageContainer (typically as a response to a
    * mouse click) and deselects the old selected container, if it
    * exists.
    * 
    * @param c the ImageContainer to select, or null to deselect
    * any currently selected container
    */
   public void setSelectedContainer(ImageContainer c) {

      if(_debug) {
	 System.out.println(">>>>>> ImageView: enter setSelectedContainer");
      }
      Set<String> sel = Collections.emptySet();
      tvim.getGeneSelection().resetCorrelationValue();

      if(c == null) {
	 if(_debug) {
	    System.out.println("passed in ImageContainer is null");
	 }
	 if (selectedImageContainer != null) {
	    selectedImageContainer.setSelected(false);
	    selectedImageContainer = c;
	    if(_debug) {
	       System.out.println("setting selectedImageContainer false");
	    }
	 }
      } else {
	 if(_debug) {
	    System.out.println("passed in ImageContainer is "+c.getID());
	 }
	 if (selectedImageContainer == null) {
	    if(_debug) {
	       System.out.println("selected ImageContainer is null");
	    }
	    selectedImageContainer = c;
	    sel = Collections.singleton(selectedImageContainer.getID());
	    if(_debug) {
	       System.out.println("setting selectedImageContainer true");
	    }
	    selectedImageContainer.setSelected(true);
	 } else {
	    if(_debug) {
	       System.out.println("selected ImageContainer is "+selectedImageContainer.getID());
	    }
	    if (c.getID().equals(selectedImageContainer.getID())) {
	       if(_debug) {
		  System.out.println("setting selectedImageContainer false");
	       }
	       selectedImageContainer.setSelected(false);
	       selectedImageContainer = null;
	    } else {
	       selectedImageContainer = c;
	       sel = Collections.singleton(selectedImageContainer.getID());
	       if(_debug) {
		  System.out.println("setting selectedImageContainer true");
	       }
	       selectedImageContainer.setSelected(true);
	    }
	 }
      }

      tvim.getGeneSelection().setHighlightedNodes(sel);
      // We only want to update the tree, so we need to prevent the
      // image view manager from being notified.
      tvim.getGeneSelection().deleteObserver(tvim);
      tvim.getGeneSelection().notifyObservers();
      tvim.getGeneSelection().addObserver(tvim);

      if(_debug) {
	 System.out.println("<<<<<< ImageView: exit setSelectedContainer");
      }
   }
//-------------------------------------------------------------------
   /** 
    * Arranges the subcomponents.
    */
   private void doMediaLayout() {
      assert top != null;
      assert bottom != null;

      //top.setVisibleRect(getVisibleRect());
      //bottom.setVisibleRect(getVisibleRect());

      if(tvim.isCorrelationView()) {
         clearTop();
      }

      setDividerSize(defaultDividerSize);
      if (top.isEmpty()) {
	 //setDividerSize(0);
	 setDividerLocation(0.0);
      } else if(bottom.isEmpty()) {
	 //setDividerSize(0);
	 setDividerLocation(1.0);
      } else {
	 //setDividerLocation(top.getPreferredSize().height);
	 setDividerLocation(0.3);
      }
      repaint();
   }
//-------------------------------------------------------------------
   /*
      public void paintComponent(Graphics g)
      {
      super.paintComponent(g);
      if (getComponentCount() == 0)
      {
      Rectangle vis = getVisibleRect();
      FontMetrics fm = g.getFontMetrics();
      String[] text = {"Click on the tree to select one branch or","move slider to select multiple branches"};
      for (int i=0; i<text.length; i++)
      g.drawString(text[i],
      vis.width/2 - fm.stringWidth(text[i])/2,
      vis.height/2 + (2*i-text.length)*fm.getHeight()/2);
      }
      }
    */

//-------------------------------------------------------------------
   /** 
    * Method from java.util.Observer that is called when an
    * observed object is changed in some way.
    *
    * The parent ImageViewManager is an observer of its
    * geneSelection object, and this class becomes involved
    * because any calls to update() in the ImageViewManager class
    * are passed to all of the TVImagePanels that it manages.
    *
    * When this happens, we need to ask the image server to load
    * some new images.
    * 
    * Most of this stuff should eventually be moved to the
    * ImageScrollPane class, as we're duplicating code here and
    * assuming that this class has exactly two subcomponents.
    * 
    * @param o 
    * @param arg 
    */
   public synchronized void update(Observable o, Object arg) {

      //_debug = true;
      if(_debug) {
	 if(o != null && arg != null) {
	    System.out.println("\n>>>>>> ImageView: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg+"------------");
	 } else if(o != null && arg == null) {
	    System.out.println("\n>>>>>> ImageView: enter update, Observable = "+o.getClass().getSimpleName()+", arg = null ------------");
	 } else if(o == null && arg != null) {
	    System.out.println("\n>>>>>> ImageView: enter update, Observable = null, arg = "+arg+"------------");
	 } else if(o == null && arg == null) {
	    System.out.println(">>>>>> ImageView: enter update, Observable = null, arg = null ......................................");
	 }
      }
      //_debug = false;

      int useCase = tvim.getUseCase();

      if (o != null && o.equals(tvim.getGeneSelection())) {

	 Collection<String> newSelectedNodes;
	 Collection<String> newSelectedGenes;
	 /*
	  * If this tab isn't capable of displaying node images but
	  * we've selected a bunch of nodes by correlation, then
	  * this tab has nothing useful to display.
	  */
	 if (!displaysSelectionByCorrelation && tvim.getGeneSelection().isCorrelationValueSet()) {
	    newSelectedNodes = Collections.emptySet();
	    newSelectedGenes = Collections.emptySet();
	 } else {
	    newSelectedNodes = tvim.getGeneSelection().getSelectedNodes();
	    newSelectedGenes = tvim.getSelectedGenes();
	 }

	 assert newSelectedNodes != null;
	 assert newSelectedGenes != null;
	 if(_debug) {
	    System.out.println("tvim.isCorrelationView() = "+tvim.isCorrelationView());
	    System.out.println("tvim.isDoubleClickGene() = "+tvim.isDoubleClickGene());
	    System.out.println("tvim.isDoubleClickNode() = "+tvim.isDoubleClickNode());
	    printCollection(newSelectedNodes, "newSelectedNodes");
	    printCollection(newSelectedGenes, "newSelectedGenes");
	 }

	 Set<String> topToAdd    = new HashSet<String>();
	 Set<String> topToDelete = new HashSet<String>();

	 Set<String> bottomToAdd    = new HashSet<String>();
	 Set<String> bottomToDelete = new HashSet<String>();

	 Set<String> newTop = new HashSet<String>();
	 Set<String> newBottom = new HashSet<String>();

	 ImageViewManager.printDebugMessage("Calculating changes in selection");

         switch(useCase) {
	    case TVTypes.CORRELATION_CHANGED:
	       newBottom.addAll(newSelectedNodes);
	       break;
	    case TVTypes.DOUBLE_CLICK_NODE:
	       newTop.addAll(newSelectedNodes);
	       break;
	    case TVTypes.DOUBLE_CLICK_GENE:
	       newTop.addAll(newSelectedGenes);
	       break;
	    case TVTypes.TREE_NODE_CLICKED:
	       newTop.addAll(newSelectedNodes);
	       break;
	    case TVTypes.MATRIX_CLICKED:
	       newTop.addAll(newSelectedGenes);
	       break;
	 } // switch
	 newBottom.addAll(newSelectedGenes);

	 for (String s : topPanels) {
	    if (!newTop.contains(s)) {
	       topToDelete.add(s);
	    }
	 }

	 for (String s : newTop) {
	    if (!topPanels.contains(s)) {
	       topToAdd.add(s);
	    }
	 }

	 topPanels.removeAll(topToDelete);
	 topPanels.addAll(topToAdd);

	 for (String s : bottomPanels) {
	    if (!newBottom.contains(s)) {
	       bottomToDelete.add(s);
	    }
	 }

	 for (String s : newBottom) {
	    bottomToAdd.add(s);
	 }

	 bottomPanels.removeAll(bottomToDelete);
	 bottomPanels.addAll(bottomToAdd);

	 top.removePanels(topToDelete);
	 bottom.removePanels(bottomToDelete);

	 setDividerSize(defaultDividerSize);
	 setDividerLocation(0.55);

	 if(_debug) {
	    printCollection(topToAdd, "topToAdd");
	    printCollection(bottomToAdd, "bottomToAdd");
	 }
	 /*
	  * If the only required operation is to remove some panels,
	  * then we have to show this on screen.  Otherwise, we'll
	  * just repaint when the new panels arrive.
	  */
	 if (topToAdd.isEmpty() && bottomToAdd.isEmpty()) {
	    doMediaLayout();
	 } else {
	    //System.out.println("\n"+viewDescription()+" is requesting new images");
	    ImageViewManager.printDebugMessage("Requesting new images");
	    if (TVTypes.MULTITHREADED) {
	       top.imgProducer.loadImagesMT(topToAdd, (Observer)this);
	       bottom.imgProducer.loadImagesMT(bottomToAdd, (Observer)this);
	    } else {
	       bottom.imgProducer.loadImagesST(bottomToAdd, this, false);
	       if(useCase != TVTypes.CORRELATION_CHANGED)  {
		  //System.out.println("\n"+viewDescription()+" is requesting new images");
		  top.imgProducer.loadImagesST(topToAdd, this, true);
	       }
	    }
	    if(bottomToAdd.isEmpty()) {
	       setDividerLocation(1.0);
	    }
	    if(topToAdd.isEmpty()) {
	       setDividerLocation(0.0);
	    }
	 }
      }
      //_debug = true;
      if(_debug) {
	 if(o != null) {
	    System.out.println("<<<<<< ImageView: exit update, Observable = "+o.getClass().getSimpleName()+"-------------------------------");
	 } else {
	    System.out.println("<<<<<< ImageView: exit update, Observable = null .....................................");
	 }
      }
      //_debug = false;
   } // update

//------------------------------------------------------------------------------------
   protected void addLoadedImages(Set<String>nodeIDs, ImageServer server, boolean isTop) {
      //System.out.println("enter addLoadedImages, isTop = "+isTop);

      if(nodeIDs == null || nodeIDs.size() <= 0) {
         System.out.println("returning because nodeIDs is empty or null");
	 return;
      }
      if(server == null) {
         System.out.println("returning because server is null");
	 return;
      }

      // top and bottom may both have the same nodeID for matrix & double click gene
      for(String str : nodeIDs) {
         if(isTop) {
	    if (topPanels.contains(str)) {
	       top.addPanel(str);
	    }
	 } else {
	    if (bottomPanels.contains(str)) {
	       bottom.addPanel(str);
	    }
	 }
      }

      doMediaLayout();

      //System.out.println("exit addLoadedImages, isTop = "+isTop);
   } // addLoadedImages

//------------------------------------------------------------------------------------
      //} else {
	 /*
	  * Assume that any other updates must be from the image
	  * server, though this may turn out to be a bad assumption.
	  *
	  * Some checking may be a good idea, but this requires
	  * identification of ImageLoaders which are currently
	  * anonymous (see ImageServer.java).
	  */
	 /*
	  * If the panel we've received isn't in the list of
	  * selected nodes (this could happen if the selection
	  * changes while images are loading), we should ignore it.
	  */
	  /*
	 boolean addedPanelToTop = false;
	 boolean addedPanelToBottom = false;
	 //System.out.println("update from image server: arg = "+(String)arg);
	 if (topPanels.contains((String)arg)) {
	    top.addPanel((String)arg);
	    //System.out.println("added to topPanels ");
	    //printCollection(topPanels,"top");
	    addedPanelToTop = true;
	 }
	 if (bottomPanels.contains((String)arg)) {
	    bottom.addPanel((String)arg);
	    //System.out.println("added to bottomPanels ");
	    //printCollection(bottomPanels,"bottom");
	    addedPanelToBottom = true;
	 }
	 //_debug = true;
	 if(!addedPanelToTop && !addedPanelToBottom) {
	    if(_debug) {
	       System.out.println("<<<<<< ImageView: update, returning no panels added .............");
	    }
	    return;
	 }

	 boolean readyToDisplay = true;

	 if (!INCREMENTAL_IMAGE_LOADING) {
	    boolean panelsContainsNode;
	    for (String s : topPanels) {
	       panelsContainsNode = false;
	       for (ImageContainer c : top.getPanels()) {
		  if (c.getID().equals(s)) {
		     panelsContainsNode = true;
		  }
	       }
	       if (!panelsContainsNode) {
		  readyToDisplay = false;
	       }
	    }
	    for (String s : bottomPanels) {
	       panelsContainsNode = false;
	       for (ImageContainer c : bottom.getPanels()) {
		  if (c.getID().equals(s)) {
		     panelsContainsNode = true;
		  }
	       }
	       if (!panelsContainsNode) {
		  readyToDisplay = false;
	       }
	    } // for
	 } else {
	     // Just to make things look smoother, we won't display
	     // anything until at least one image in the top panel has
	     // been loaded
	    if (top.getPanels().isEmpty() && !topPanels.isEmpty()) {
	       readyToDisplay = false;
	    }
	 }

	 if (readyToDisplay) {
	    doMediaLayout();
	 }
      } // else
      */

   private void printCollection(Collection col, String name) {
      if(col == null || col.size() <= 0) {
	 System.out.println("Collection "+name+" ----------------- contains 0 entries");
	 return;
      } else {
	 System.out.println("Collection "+name+" ----------------- contains "+col.size()+" entries:");
      }

      Iterator it = col.iterator();
      while(it.hasNext()) {
         System.out.println(it.next().toString());
      }
   }

   /**
    *   We want to remove any gene images that may have been put in the top panel if it is a 'correlationView'
    */
   protected void clearTop() {
      if(topPanels == null || topPanels.size() <= 0) {
         return;
      }
      top.getContent().removeAll();
      topPanels.clear();
   }

   /**
    *   Clear all ImageContainers from memory.
    */
   protected void clearAll() {
      if(topPanels != null && topPanels.size() > 0) {
	 topPanels.clear();
      }
      if(bottomPanels != null && bottomPanels.size() > 0) {
	 bottomPanels.clear();
      }
      top.clearAll();
      bottom.clearAll();
   }

} // class ImageView
