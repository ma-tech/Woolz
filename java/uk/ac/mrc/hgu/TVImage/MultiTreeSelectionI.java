package uk.ac.mrc.hgu.TVImage;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;
import java.util.*;

/**
 * An extension of TreeSelectionI that allows multiple nodes to be
 * selected.
 *
 * A single node selection is represented by a single element in
 * the ArrayList of selected nodes.  If multiple nodes are selected,
 * getSelectedNode() should return null.
 *
 * The intention is that existing TreeSelection objects may simply
 * be cast to a MultiTreeSelection where required.
 */
public interface MultiTreeSelectionI extends TreeSelectionI
{
	public void printDetails(String requester);
	public void setLeftTreeDrawer(LeftTreeDrawer d);

	/**
	 * Selects tree nodes.
	 *
	 * @param n IDs of nodes to select
	 */
	public void setSelectedNodes(Set<String> n);

	public void setNode(String n, boolean selected);
	public boolean isSelected(String n);
	
	/**
	 * Gets the selected tree nodes
	 *
	 * @return nodeIDs of selected nodes 
	 */
	public Set<String> getSelectedNodes();

	/** 
	 * Returns the current correlation value, used to select multiple
	 * nodes in the tree according to their correlation.
	 * 
	 * @return 
	 */
	public double getCorrelationValue();
	
	/** 
	 * Set the correlation value, as described above.
	 * 
	 * @param d 
	 */
	public void setCorrelationValue(double d);

	public void resetCorrelationValue();
	public boolean isCorrelationValueSet();

	public void setHighlightedNodes(Set<String> n);
	public Set<String> getHighlightedNodes();

	public abstract void notifyObservers(Object arg);
}
