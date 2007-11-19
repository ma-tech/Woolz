package uk.ac.mrc.hgu.TVImage;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;
import java.util.*;

/**
 * An extension of TreeSelection that allows multiple nodes to be
 * selected.
 *
 * Further details available in the interface documentation
 * (MultiTreeSelectionI.java)
 *
 * @see uk.ac.mrc.hgu.TVImage.MultiTreeSelectionI
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 * @version
 */
public class MultiTreeSelection extends TreeSelection implements MultiTreeSelectionI {

   // What we really want is something that's both a list and a
   // set (not a SortedSet), as it should preserve ordering based
   // on the order in which elements are added, rather than by
   // some specified comparator.
   //
   // For now, we'll just use a standard set.
   private Set<String> selectedNodes = new HashSet<String>();
   private Set<String> highlightedNodes = new HashSet<String>();
   private double correlation = -1;
   private boolean _debug = false;

   // XXX: Total hack which allows us to examine the structure of
   // the tree.  We don't use it to draw anything.
   // Perhaps we can get this information from GTRInfo instead?
   /*very*/ private LeftTreeDrawer drawer = null;

   private void synchMap() {
      if(_debug) {
	 System.out.println(">>>>>> MultiTreeSelection: enter synchMap");
      }
      if (selectedNodes != null) {
	 double corr = getCorrelationValue();
	 //System.out.println("correlation = "+corr);
	 if (corr < 0) {
	    for (String n : selectedNodes) {
	       TreeDrawerNode tdn = drawer.getNodeById(n);
	       int start = (int)tdn.getLeftLeaf().getIndex();
	       int end   = (int)tdn.getRightLeaf().getIndex();
	       selectIndexRange(start, end);
	    }
	 } else {
	    //System.out.println("selecting all indexes");
	    selectAllIndexes();

	    for (String n : selectedNodes) {
	       TreeDrawerNode tdn = drawer.getNodeById(n);
	       int start = (int)tdn.getLeftLeaf().getIndex();
	       int end   = (int)tdn.getRightLeaf().getIndex();
	       for (int i=start; i<=end; i++) {
		  setIndex(i,false);
	       }
	    }
	 }
      }
      if(_debug) {
	 System.out.println("<<<<<< MultiTreeSelection: exit synchMap");
      }
   }

   /**
    * Constructor for the MultiTreeSelection object
    *
    * @param nIndex number of indexes which can be selected
    */
   public MultiTreeSelection(int nIndex) {
      super(nIndex);
   }

   public void setLeftTreeDrawer(LeftTreeDrawer d) {
      drawer = d;
   }

   public void printDetails(String requester) {
      System.out.println("printDetails() called from "+requester);
      Iterator<String> it;
      it = selectedNodes.iterator();
      while (it.hasNext()) {
	 System.out.println("Node: " + it.next());
      }
      it = highlightedNodes.iterator();
      while (it.hasNext()) {
	 System.out.println("H-Node: " + it.next());
      }
      int[] indexes = getSelectedIndexes();
      for (int i=0; i<indexes.length; i++) {
	 System.out.println("Index " + i + ": " + indexes[i]);
      }
      System.out.println("MinIndex: " + getMinIndex());
      System.out.println("MaxIndex: " + getMaxIndex());
   }

   /* (non-Javadoc)
    * @see edu.stanford.genetics.treeview.TreeSelectionI#setSelectedNode(java.lang.String)
    */
   public void setSelectedNode(String n) {
      if(_debug) {
	 System.out.println(">>>>>> MultiTreeSelection: enter setSelectedNode");
	 System.out.println(n);
      }
      Set<String> s;
      if (n == null) {
	 s = Collections.emptySet();
      } else {
	 s = Collections.singleton(n);
      }
      setSelectedNodes(s);
      if(_debug) {
	 System.out.println(">>>>>> MultiTreeSelection: exit setSelectedNode");
      }
   }

   /* (non-Javadoc)
    * @see uk.ac.mrc.hgu.TVImage.MultiTreeSelectionI#setSelectedNodes()
    */
   public void setSelectedNodes(Set<String> n) {
      if(_debug) {
	 System.out.println(">>>>>> MultiTreeSelection: enter setSelectedNodes");
	 if(n == null) {
	    System.out.println("n is null");
	 } else if(n.size() <=0) {
	    System.out.println("n is empty");
	 } else {
	    for(String s : n) {
	       System.out.println(s);
	    }
	 }
      }
      //System.out.println("selectedNodes.clear()");
      selectedNodes.clear();
      if (n != null) {
	 //System.out.println("selectedNodes.addAll(n)");
	 selectedNodes.addAll(n);
      }
      setHighlightedNodes(n);
      ImageViewManager.printDebugMessage("Clearing selected leaf list");
      deselectAllIndexes();
      ImageViewManager.printDebugMessage("Selecting constituent leaves for selected node");
      synchMap();
      //System.out.println("setSelectedNodes: about to setChanged() **************************");
      setChanged();
      if(_debug) {
	 System.out.println("<<<<<< MultiTreeSelection: exit setSelectedNodes");
      }
   }

   public void setNode(String n, boolean s) {
      if(_debug) {
	 if ((s && selectedNodes.contains(n)) || (!s && !selectedNodes.contains(n))) {
	    return;
	 } else if (s) {
	    selectedNodes.add(n);
	 } else {
	    selectedNodes.remove(n);
	 }
      }
      setHighlightedNodes(selectedNodes);
      deselectAllIndexes();
      synchMap();
      //System.out.println("setNode: about to setChanged() **************************");
      setChanged();
      if(_debug) {
	 System.out.println("<<<<<< MultiTreeSelection: exit setSelectedNodes");
      }
   }

   public boolean isSelected(String n) {
      return (selectedNodes.contains(n));
   }

   /* (non-Javadoc)
    * @see edu.stanford.genetics.treeview.TreeSelectionI#getSelectedNode()
    */
   public String getSelectedNode() {
      Set<String> nodes = getSelectedNodes();
      if (nodes.size() == 1) {
	 return (nodes.toArray(new String[0])[0]);
      } else {
	 return null;
      }
   }

   /* (non-Javadoc)
    * @see uk.ac.mrc.hgu.TVImage.MultiTreeSelectionI#getSelectedNodes()
    */
   public Set<String> getSelectedNodes() {
      return selectedNodes;
   }

   public double getCorrelationValue() {
      return correlation;
   }

   public void setCorrelationValue(double d) {
      if(_debug) {
	 System.out.println(">>>>>> MultiTreeSelection: enter setCorrelationValue "+d);
      }
      if (d == correlation) {
	 if(_debug) {
	    System.out.println("<<<<<< MultiTreeSelection: returning, d == correlation "+d);
	 }
	 return;
      }
      correlation = d;
      //System.out.println("setCorrelationValue: about to setChanged() **************************");
      setChanged();
      if(_debug) {
	 System.out.println("<<<<<< MultiTreeSelection: exit setCorrelationValue "+d);
      }
   }

   public void resetCorrelationValue() {
      setCorrelationValue(-1);
   }

   public boolean isCorrelationValueSet() {
      return (getCorrelationValue() != -1);
   }

   public void setHighlightedNodes(Set<String> n) {
      if(_debug) {
	 System.out.println(">>>>>> MultiTreeSelection: enter setHighlightedNodes");
	 if(n == null) {
	    System.out.println("set is null");
	 } else if(n.size() <= 0) {
	    System.out.println("set is empty");
	 } else {
	    for(String s : n) {
	       System.out.println(s);
	    }
	 }
      }
      if (n != null) {
	 highlightedNodes = n;
      } else {
	 highlightedNodes.clear();
      }
      //System.out.println("setHighlightedNodes: about to setChanged() **************************");
      setChanged();
      if(_debug) {
	 System.out.println("<<<<<< MultiTreeSelection: exit setHighlightedNodes");
      }
   }

   public Set<String> getHighlightedNodes() {
      return highlightedNodes;
   }
}
