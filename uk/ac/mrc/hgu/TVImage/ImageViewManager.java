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
 * Manager for image-related views.  Also contains a bunch of
 * wrapper methods for the geneSelection and nodeSelection
 * objects.
 *
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 * @version 
 */
public class ImageViewManager extends JTabbedPane implements Observer, TVTypes {

   private MultiTreeSelectionI geneSelection = null;
   private JTabbedPane tabbedPane = null;
   private HeaderInfo gtrInfo = null;
   private HeaderInfo geneInfo = null;
   private ArrayList<JComponent> views = new ArrayList<JComponent>();
   private ArrayList<String> viewNames = new ArrayList<String>();
   private int scrollValue = 0;
   private int useCase = NONE;
   private XmlAppSettings appSettings;
   private ScrollGroup sg = new ScrollGroup();
   private ModelManipulator mm;
   private ViewFrame viewFrame;
   private LeftTreeDrawer drawer;
   private TreeDrawerNode rootNode;
   private Object _updateParam = null;
   private boolean isCorrelationView = false;
   private boolean isDoubleClickGene = false;
   private boolean isDoubleClickNode = false;
   private boolean tabChangedByUser = true;

   private boolean _debug = false;
   private static final boolean DEBUG_MODE = false;

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

      this.addChangeListener(new ChangeListener() {
	 public void stateChanged(ChangeEvent e) {
	    doTabChanged(e);
	 }
      });

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
      ImageScrollPane nodePane
	 = new ImageScrollPane(new ScrollGroup(), NODE, hif, this);
      ImageScrollPane nodePane2
	 = new ImageScrollPane(new ScrollGroup(), NODE, hif, this);
      ImageScrollPane nodePane3
	 = new ImageScrollPane(new ScrollGroup(), NODE, hif, this);
      ImageScrollPane geneHeatmapPane
	 = new ImageScrollPane(sg, HEATMAP,   hif, this);
      ImageScrollPane geneRawPane
	 = new ImageScrollPane(sg, RAW,       rif, this);
      ImageScrollPane geneAnnotatedPane
	 = new ImageScrollPane(sg, ANNOTATED, aif, this);

      ImageView heatmapView = new ImageView(
	    this,
	    appSettings,
	    true,
	    HEATMAP_DESCRIPTION,
	    nodePane,
	    geneHeatmapPane);
      ImageView rawView = new ImageView(
	    this,
	    appSettings,
	    false,
	    RAW_DESCRIPTION,
	    nodePane2,
	    geneRawPane);
      ImageView annotatedView = new ImageView(
	    this,
	    appSettings,
	    false,
	    ANNOTATED_DESCRIPTION,
	    nodePane3,
	    geneAnnotatedPane);

      views.add(heatmapView);
      views.add(rawView);
      views.add(annotatedView);
      views.add(new SummaryView(this));

      viewNames.add(heatmapView.viewDescription());
      viewNames.add(rawView.viewDescription());
      viewNames.add(annotatedView.viewDescription());
      viewNames.add(SUMMARY_DESCRIPTION);

      for (int i=0; i<views.size(); i++) {
	 if(_debug) {
	    System.out.println("ImageViewManager: adding "+viewNames.get(i));
	 }
	 addTab(viewNames.get(i), views.get(i));
      }
   } // constructor

   public int getUseCase() {
      return useCase;
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
      _updateParam = arg;
      if (o == getGeneSelection()) {
	 if (getGeneSelection().isCorrelationValueSet()) {
	    /*
	     * XXX: we're currently hard coding allowable tab indexes
	     * for correlation views.  Bad.
	     */
	    tabChangedByUser = false;
	    setSelectedIndex(0);
	    tabChangedByUser = true;
	    setEnabledAt(1, false);
	    setEnabledAt(2, false);
	    setEnabledAt(3, false);
	    setCorrelationView(true);
	 } else {
	    tabChangedByUser = false;
	    setSelectedIndex(1);
	    tabChangedByUser = true;
	    setEnabledAt(1, true);
	    setEnabledAt(2, true);
	    setEnabledAt(3, true);
	    setCorrelationView(false);
	 }
	 clearAll();
	  // We only want to update views that will be displayed.
	 if (isCorrelationView()) {
	    useCase = CORRELATION_CHANGED;
	    updateView(o, arg, HEATMAP_DESCRIPTION);
	 } else if (arg != null && arg.toString().trim().equalsIgnoreCase("treeNodeClicked")) {
	    useCase = TREE_NODE_CLICKED;
	    updateView(o, arg, RAW_DESCRIPTION);
	    updateView(o, arg, SUMMARY_DESCRIPTION);
	 } else if (isDoubleClickNode()) {
	    useCase = DOUBLE_CLICK_NODE;
	    updateView(o, arg, RAW_DESCRIPTION);
	    updateView(o, arg, SUMMARY_DESCRIPTION);
	 } else if (isDoubleClickGene()) {
	    useCase = DOUBLE_CLICK_GENE;
	    updateView(o, arg, RAW_DESCRIPTION);
	    updateView(o, arg, SUMMARY_DESCRIPTION);
	 } else {
	    useCase = MATRIX_CLICKED;
	    updateView(o, arg, RAW_DESCRIPTION);
	    updateView(o, arg, SUMMARY_DESCRIPTION);
	 }
      } else if (o == drawer && (rootNode == null || !rootNode.equals(drawer.getRootNode()))) {
	 rootNode = drawer.getRootNode();
	 mm.doTreeSearch(rootNode);
      }
      if(_debug) {
	 System.out.println("<<<<<< ImageViewManager: exit update, Observable = "+o.getClass().getSimpleName());
      }
   } // update

   /**
    *   Update the relevant image view.
    *   @param o the Observable object that caused this update.
    *   @param arg Data to be passed to update.
    *   @param description the type ov view to be updated.
    */
   public void updateView(Observable o, Object arg, String description) {
      if(views == null || views.size() <= 0) {
         return;
      }
      int siz = views.size();
      JComponent view;
      for (int i=0; i<siz; i++) {
	 view = views.get(i);
	 if(getViewDescription(view).equalsIgnoreCase(description)) {
	    //System.out.println("ImageViewManager: updating "+description);
	    ((Observer)view).update(o,arg);
	    break;
	 }
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

   public void clearAll() {
      //System.out.println("ImageViewManager: clearAll");
      if(views == null) {
         return;
      }
      int siz = views.size();
      if(siz <= 0) {
         return;
      }
      JComponent view;
      // SummaryView is a different class from the other views.
      for (int i=0; i<siz; i++) {
	 view = views.get(i);
	 if(getViewDescription(view).equalsIgnoreCase(TVTypes.SUMMARY_DESCRIPTION)) {
	    continue;
	 }
	 ((ImageView)views.get(i)).clearAll();
      }
   } // clearAll

   //---------------------------------------------------------------------------
   public static String getViewDescription(Object view) {

      String ret = "";
      Method M1 = null;

      try {
	 M1 = view.getClass().getMethod("viewDescription", null);
	 ret = (String)M1.invoke(view, null);
      }
      catch (InvocationTargetException ie) {}
      catch (NoSuchMethodException ie) {}
      catch (IllegalAccessException ie) {}

      return ret;
   } // getViewDescription

   //---------------------------------------------------------------------------
   private void doTabChanged(ChangeEvent e) {
      int index = getSelectedIndex();
      String title = getTitleAt(index);
      if(!tabChangedByUser) {
	 return;
      }
      clearAll();
      Observable o = (Observable)getGeneSelection();
      JComponent view;
      for (int i=0; i<views.size(); i++) {
	 view = views.get(i);
	 if(getViewDescription(view).equalsIgnoreCase(title)) {
	    ((Observer)view).update(o, _updateParam);
	 }
      }
   } // doTabChanged

}
