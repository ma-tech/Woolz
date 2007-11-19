package uk.ac.mrc.hgu.TVImage;

import java.util.*;
import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/** 
 * Some modifications that we wish to perform on the
 * geneSelection object can only be done with some details held
 * in geneInfo, nodeInfo etc.  Here, we provide some wrappers
 * that allow us to make these modifications.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 * @version 
 */
class ModelManipulator
{
	private MultiTreeSelectionI geneSelection = null;
	private HeaderInfo nodeInfo = null;
	private HeaderInfo geneInfo = null;
	private HeaderSummary headerSummary = new HeaderSummary();
	private List<String> inorderTreeSearch = new LinkedList<String>();

	private boolean _debug = false;

	public int getInorderSearchIndex(String node)
	{ return inorderTreeSearch.indexOf(node); }

	public int getInorderSearchSize()
	{ return inorderTreeSearch.size(); }

	public ModelManipulator(HeaderInfo nI, HeaderInfo gI)
	{
		nodeInfo = nI;
		geneInfo = gI;
	}

	public void doTreeSearch(TreeDrawerNode root)
	{
		//inorderTreeSearch.clear();
		inorderTreeSearch = new LinkedList<String>();
		Stack<TreeDrawerNode> remaining = new Stack<TreeDrawerNode>();
		TreeDrawerNode temp = root;

		do
		{
			while (temp != null)
			{
				remaining.push(temp);
				temp = temp.getLeft();
			}

			if (!remaining.empty())
				temp = (TreeDrawerNode)remaining.pop();

			if (temp != null)
			{
				inorderTreeSearch.add(temp.getId());
				temp = temp.getRight();
			}
		} while (!remaining.empty() || temp != null);
	}

	public String getEMAGEID(int index)
	{
		String out = null;
		String[] summaryArray = getSummaryArray(index);
		if ((summaryArray == null))
			return null;
		return summaryArray[0];
	}

	/**
	 * test to see what is stored for a gene header
	 *
	 */
	public String getEMAGEID(String id)
	{
		int index = getGeneIndex(id);
		return getEMAGEID(index);
	}

	/**
	 * Display GENEX ID in GENEX:123 format
	 *
	 * @param index 
	 * @return 
	 */
	public String getGeneName(int index)
	{
		if (index >= 0)
			return geneInfo.getHeader(index, "NAME");
		else
			return "";
	}

	/**
	 * test to see what is stored for a gene header
	 *
	 */
	public String getGeneName(String id)
	{
		int index = getGeneIndex(id);
		return getGeneName(index);
	}


	public String[] getSummaryArray(int index)
	{
		if (geneInfo == null)
			return null;
		if (index < 0 || index > geneInfo.getNumHeaders()-1)
			return null;
		return headerSummary.getSummaryArray(geneInfo, index);
	}

	public int[] getSelectedIndexes()
	{
		return geneSelection.getSelectedIndexes();
	}

	public Set<String> getSelectedGenes() {
	   int[] selectedIndexes = getSelectedIndexes();
	   Set<String> selectedGenes = new HashSet<String>();
	   for (int i=0; i<selectedIndexes.length; i++) {
	      //System.out.println("getSelectedGenes: selectedIndexes["+i+"] = "+selectedIndexes[i]);
	      if (selectedIndexes[i] >= 0) {
		 selectedGenes.add(getNodeID(selectedIndexes[i]));
		 //System.out.println("getSelectedGenes: adding "+getNodeID(selectedIndexes[i]));
	      }
	   }
	   return selectedGenes;
	}

	public String getNodeID(int index)
	{ return geneInfo.getHeader(index,"GID"); }

	public int getGeneIndex(String nodeID)
	{ return geneInfo.getHeaderIndex(nodeID); }

	public int getNodeIndex(String nodeID)
	{ return nodeInfo.getHeaderIndex(nodeID); }

	public int getMinGeneIndex() {
	   if (geneSelection != null) {
	      return geneSelection.getMinIndex();
	   } else {
	      return 0;
	   }
	}

	public String getCorrelation(String nodeID) {
	   if (nodeInfo == null || nodeID == null) {
	      return null;
	   }
	   int headerIndex = getNodeIndex(nodeID);
	   String correlation = null;
	   if (headerIndex >= 0) {
	      correlation = nodeInfo.getHeader(headerIndex, "CORRELATION");
	   }
	   return correlation;
	}

	public void selectNodesByCorrelation(double corr) {
	   if(_debug) {
	      System.out.println(">>>>>> ModelManipulator: enter selectNodesByCorrelation "+corr);
	   }
	   geneSelection.setCorrelationValue(corr);
	   // Root node is always(?) the last node in the list
	   int index = nodeInfo.getNumHeaders() - 1;
	   String rootNode = nodeInfo.getHeader(index, "NODEID");

	   if(_debug) {
	      System.out.println("index = "+index);
	      System.out.println("rootNode = "+rootNode);
	   }

	   Set<String> correlatedNodes = new HashSet<String>();
	   Stack<String> remaining = new Stack<String>();
	   remaining.push(rootNode);
	   //System.out.println("enter while loop");
	   while (remaining.empty() == false) {
	      String node = remaining.pop();
	      assert node != null;
	      //System.out.println("popped "+node+" from remaining");

	      index = getNodeIndex(node);
	      //System.out.println(node+"'s index = "+index);

	      // If not a leaf...
	      if (index >= 0) {
		 String correlation = getCorrelation(node);
		 assert correlation != null;
		 //System.out.println("   correlation = "+correlation);

		 if (Double.parseDouble(correlation) > corr) {
		    correlatedNodes.add(node);
		    //System.out.println("add "+node+" to correlatedNodes");
		 } else {
		    String left = nodeInfo.getHeader(index, "LEFT");
		    String right = nodeInfo.getHeader(index, "RIGHT");
		    remaining.push(left);
		    remaining.push(right);
		    //System.out.println("pushed "+left+" and "+right+" onto remaining");
		 }
	      }
	   }
	   //System.out.println("exit while loop");

	   geneSelection.setSelectedNodes(correlatedNodes);
	   geneSelection.notifyObservers();
	   if(_debug) {
	      System.out.println(">>>>>> ModelManipulator: exit selectNodesByCorrelation "+corr);
	   }
	}

	public MultiTreeSelectionI getGeneSelection() {
	   return geneSelection;
	}

	public void setGeneSelection(TreeSelectionI geneSelection) {
	   this.geneSelection = (MultiTreeSelectionI)geneSelection;
	}
}
