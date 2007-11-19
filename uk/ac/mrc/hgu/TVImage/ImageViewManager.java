/* file name  : ImageViewManager.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>, Paul Smith <psmith@hgu.mrc.ac.uk>
 * created    : Fri 13 Jul 2007 16:05:11 BST
 * copyright  : 
 *
 * modifications:
 * Paul Smith - fixed functionality to find gene names
 */

package uk.ac.mrc.hgu.TVImage;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.event.*;
import java.util.*;
import javax.imageio.*;
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.lang.reflect.*;
import java.net.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/** 
 * Manager for image-related panels.  Also contains a bunch of
 * wrapper methods for the geneSelection and nodeSelection
 * objects.
 *
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 * @version 
 */
public class ImageViewManager extends JTabbedPane implements Observer {

   private MultiTreeSelectionI geneSelection = null;
   private JTabbedPane tabbedPane = null;
   private HeaderInfo gtrInfo = null;
   private HeaderInfo geneInfo = null;
   private ArrayList<JComponent> panels = new ArrayList<JComponent>();
   private ArrayList<String> panelNames = new ArrayList<String>();
   private int scrollValue = 0;
   private XmlAppSettings appSettings;
   private ScrollGroup sg = new ScrollGroup();
   private ModelManipulator mm;
   private ViewFrame viewFrame;
   private LeftTreeDrawer drawer;
   private TreeDrawerNode rootNode;
   private boolean isCorrelationView = false;
   private boolean isDoubleClickGene = false;
   private boolean isDoubleClickNode = false;
   private boolean _debug = false;

   public static final Color BGCOLOR = Color.white;
   public static final boolean MULTITHREADED = false;
   public static final boolean PRELOAD = false;
   private static final boolean DEBUG_MODE = false;

   /*
    * These are bad, but they ensure that only one scroll pane of
    * a particular type can be included in a scroll group.
    */
   public static final int TYPE_NULL = 0;
   public static final int TYPE_NODE = 1;
   public static final int TYPE_HEATMAP = 2;
   public static final int TYPE_RAW = 3;
   public static final int TYPE_ANNOTATED = 4;

   public static boolean isDebugMode() {
      return DEBUG_MODE;
   }

   public static void printDebugMessage(String s) {
      if (DEBUG_MODE) System.out.println(s);
   }

   /** 
    * Construct a new ImageViewManager.  Subpanels are created and
    * added to the list of available tabs.
    */
   public ImageViewManager(HeaderInfo gtrI, HeaderInfo geneI,
	 String xmlLoc, ViewFrame vf, LeftTreeDrawer d) {

      super();

      if(_debug) {
	 System.out.println("ImageViewManager: xmlLoc = "+xmlLoc);
      }

      try {
	 appSettings = new XmlAppSettings(xmlLoc, "TreeView");
      }
      catch (Exception e) {
	 System.err.println("Could not parse XML file");
	 return;
      }

      gtrInfo = gtrI;
      geneInfo = geneI;
      viewFrame = vf;

      mm = new ModelManipulator(gtrI, geneI);
      drawer = d;
      drawer.addObserver(this);

      HeatmapImageFactory hif
	 = new HeatmapImageFactory  (this, appSettings);
      RawImageFactory rif
	 = new RawImageFactory      (this, appSettings);
      AnnotatedImageFactory aif
	 = new AnnotatedImageFactory(this, appSettings);

      /*
       * It seems that we can't add the same scroll pane to many
       * image views, so we have to create a few identical scroll
       * panes.  They use the same image factory (and image cache)
       * though, so it's not so bad.
       */
      TopImageScrollPane nodePane
	 = new TopImageScrollPane(new ScrollGroup(), TYPE_NODE, hif, this);
      TopImageScrollPane nodePane2
	 = new TopImageScrollPane(new ScrollGroup(), TYPE_NODE, hif, this);
      TopImageScrollPane nodePane3
	 = new TopImageScrollPane(new ScrollGroup(), TYPE_NODE, hif, this);
      ImageScrollPane geneHeatmapPane
	 = new ImageScrollPane(sg, TYPE_HEATMAP,   hif, this);
      ImageScrollPane geneRawPane
	 = new ImageScrollPane(sg, TYPE_RAW,       rif, this);
      ImageScrollPane geneAnnotatedPane
	 = new ImageScrollPane(sg, TYPE_ANNOTATED, aif, this);

      ImageView heatmapView = new ImageView(
	    this,
	    appSettings,
	    true,
	    "Heat map",
	    nodePane,
	    geneHeatmapPane);
      ImageView rawView = new ImageView(
	    this,
	    appSettings,
	    false,
	    "Raw data",
	    nodePane2,
	    geneRawPane);
      ImageView annotatedView = new ImageView(
	    this,
	    appSettings,
	    false,
	    "Annotated images",
	    nodePane3,
	    geneAnnotatedPane);

      panels.add(heatmapView);
      panels.add(rawView);
      panels.add(annotatedView);
      panels.add(new SummaryView(this));

      panelNames.add(heatmapView.viewDescription());
      panelNames.add(rawView.viewDescription());
      panelNames.add(annotatedView.viewDescription());
      panelNames.add(new String("Summary"));

      for (int i=0; i<panels.size(); i++) {
	 if(_debug) {
	    System.out.println("ImageViewManager: adding "+panelNames.get(i));
	 }
	 addTab(panelNames.get(i), panels.get(i));
      }
   }

   /*
    * Yuck.  Fix these.
    */
   public void displayURL(String url) {
      viewFrame.displayURL(url); 
   }

   public void displayURL(URL url) {
      viewFrame.displayURL(url.toString()); 
   }

   public String getEMAGEID(int index) {
      return mm.getEMAGEID(index); 
   }

   public String getEMAGEID(String id) {
      return mm.getEMAGEID(id); 
   }

   public String getGeneName(int index) {
      return mm.getGeneName(index); 
   }

   public String getGeneName(String id) {
      return mm.getGeneName(id); 
   }

   public String[] getSummaryArray(int index) {
      return mm.getSummaryArray(index); 
   }

   public int[] getSelectedIndexes() {
      return mm.getSelectedIndexes(); 
   }

   public Set<String> getSelectedGenes() {
      return mm.getSelectedGenes(); 
   }

   public String getNodeID(int index) {
      return mm.getNodeID(index); 
   }

   public int getGeneIndex(String nodeID) {
      return mm.getGeneIndex(nodeID); 
   }

   public int getNodeIndex(String nodeID) {
      return mm.getNodeIndex(nodeID); 
   }

   public int getMinGeneIndex() {
      return mm.getMinGeneIndex(); 
   }

   public String getCorrelation(String nodeID) {
      return mm.getCorrelation(nodeID); 
   }

   public void selectNodesByCorrelation(double corr) {
      mm.selectNodesByCorrelation(corr); 
   }

   public int getInorderSearchIndex(String node) {
      return mm.getInorderSearchIndex(node); 
   }

   public int getNumNodes() {
      //return mm.getInorderSearchSize();
      return gtrInfo.getNumHeaders();
   }

   /**
    * From java.util.Observer.
    *
    * @param o @param arg 
    */
   public void update(Observable o, Object arg) {
      if(_debug) {
	 System.out.println(">>>>>> ImageViewManager: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
      }
      if (o == getGeneSelection()) {
	 if (getGeneSelection().isCorrelationValueSet()) {
	    /*
	     * XXX: we're currently hard coding allowable tab indexes
	     * for correlation views.  Bad.
	     */
	    setSelectedIndex(0);
	    setEnabledAt(1, false);
	    setEnabledAt(2, false);
	    setEnabledAt(3, false);
	    setCorrelationView(true);
	 } else {
	    setSelectedIndex(1);
	    setEnabledAt(1, true);
	    setEnabledAt(2, true);
	    setEnabledAt(3, true);
	    setCorrelationView(false);
	 }
      } else if (o == drawer && (rootNode == null || !rootNode.equals(drawer.getRootNode()))) {
	 rootNode = drawer.getRootNode();
	 mm.doTreeSearch(rootNode);
      }

      /*
       * We may as well tell everyone what's going on (unless
       * overhead becomes a problem) -- chain of responsibility
       * design pattern?
       *
       * Also a good idea because we shouldn't expect the manager
       * instance to know details of its subpanels.
       */
      for (int i=0; i<panels.size(); i++) {
	 ((Observer)panels.get(i)).update(o,arg);
      }
      if(_debug) {
	 System.out.println("<<<<<< ImageViewManager: exit update, Observable = "+o.getClass().getSimpleName());
      }
   }

   public MultiTreeSelectionI getGeneSelection() {
      return mm.getGeneSelection();
   }

   public void setGeneSelection(TreeSelectionI geneSelection) {
      if(_debug) {
	 System.out.println("ImageViewManager: enter setGeneSelection");
      }
      if (this.geneSelection != null) {
	 this.geneSelection.deleteObserver(this);
      }
      this.geneSelection = (MultiTreeSelectionI)geneSelection;
      //this.geneSelection.printDetails("ImageViewManager.setGeneSelection");
      mm.setGeneSelection(geneSelection);
      if (this.geneSelection != null) {
	 this.geneSelection.addObserver(this);
      }
      if(_debug) {
	 System.out.println("ImageViewManager: exit setGeneSelection");
      }
   }


   public XmlAppSettings getAppSettingsObject() {
      return appSettings;
   }

   /*
    *   true if the correlation box or slider has changed and we haven't subsequently
    *   double clicked an image, or selected a tree node, or clicked in the GlobalView.
    */
   public boolean isCorrelationView() {
      return isCorrelationView;
   }

   public void setCorrelationView(boolean bool) {
      if(_debug) {
	 System.out.println("IVM: setCorrelationView "+bool);
      }
      isCorrelationView = bool;
   }

   /*
    *   true if we have double clicked a node image
    */
   public boolean isDoubleClickNode() {
      return isDoubleClickNode;
   }

   /*
    *   true if we have double clicked a gene image
    */
   public boolean isDoubleClickGene() {
      return isDoubleClickGene;
   }

   public void setDoubleClickNode(boolean bool) {
      isDoubleClickNode = bool;
   }

   public void setDoubleClickGene(boolean bool) {
      isDoubleClickGene = bool;
   }

   public void setDoubleClick(boolean bool) {
      isDoubleClickNode = bool;
      isDoubleClickGene = bool;
   }

}
