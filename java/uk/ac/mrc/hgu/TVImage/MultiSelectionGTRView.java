package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/** 
 * Modifies GTRView functionality such that multiple nodes can be selected.
 * 
 * @author Tom Perry
 * @version 
 */
public class MultiSelectionGTRView extends GTRView implements MouseListener, KeyListener {

	protected MultiTreeSelectionI geneSelection = null;
	//protected TreeDrawerNode[] currentSelectedNodes = null;
	//protected TreeDrawerNode[] currentHighlightedNodes = null;
	protected MultiSelectionLeftTreeDrawer drawer = null;
	private boolean _debug = false;

	public MultiSelectionGTRView()
	{
		super();
	}

	/** 
	 * Set the drawer
	 *
	 * @param d    The new drawer
	 */
	public void setLeftTreeDrawer(LeftTreeDrawer d)
	{
		if (drawer != null)
			drawer.deleteObserver(this);
		drawer = (MultiSelectionLeftTreeDrawer)d;
		drawer.addObserver(this);
	}

	/**
	 * Set geneSelection
	 *
	 * @param geneSelection The TreeSelection which is set by selecting
	 * genes in the GlobalView
	 */
	public void setGeneSelection(TreeSelectionI geneSelection)
	{
		if (this.geneSelection != null)
			this.geneSelection.deleteObserver(this);
		this.geneSelection = (MultiTreeSelectionI)geneSelection;
		this.geneSelection.addObserver(this);
	}

	/** 
	 * Synchronizes selected indexes with selected nodes.
	 *
	 * Ideally, this method would be in a class more directly tied
	 * to the model, but the only way to convert TreeDrawerNodes to
	 * indexes appears to be through a TreeDrawer object, which is
	 * provided here.
	 *
	 * Notifying observers of the geneSelection once to signal a
	 * node change and once more to signal a gene change is
	 * seriously not cool.
	 * 
	 * @param ns 
	 */
	 /*
	 * (nickb: it doesn't appear to be called by anything)
	private void synchMap(TreeDrawerNode[] ns) {

	   System.out.println(">>>>>> MultiSelectionGTRView: enter synchMap");

	   // According to new update spec, geneSelection has to have
	   // actually changed for this to be called.

	   if (ns != null) {
	      double corr = geneSelection.getCorrelationValue();
	      System.out.println("Correlation = "+corr);
	      if (corr < 0) {
		 for (int i=0; i<ns.length; i++) {
		    int start = (int)ns[i].getLeftLeaf().getIndex();
		    int end   = (int)ns[i].getRightLeaf().getIndex();
		    geneSelection.selectIndexRange(start, end);
		 }
	      } else {
		 geneSelection.selectAllIndexes();

		 for (TreeDrawerNode n : ns) {
		    int start = (int)n.getLeftLeaf().getIndex();
		    int end   = (int)n.getRightLeaf().getIndex();
		    for (int i=start; i<=end; i++)
		       geneSelection.setIndex(i,false);
		 }
	      }
	   }
	   System.out.println("<<<<<< MultiSelectionGTRView: exit synchMap");
	}
	 */

	public void redrawTree(TreeDrawerNode[] ns)
	{
		double corr = geneSelection.getCorrelationValue();
		if (getYScaleEq() != null)
		{
			/*
			 * Selective repainting is no longer possible because increasing
			 * the correlation leaves white dots where the correlation line
			 * used to be
			 *
			//Repaint all subtrees that are no longer selected
			if (ons != null)
				for (int i=0; i<ons.length; i++)
					if (!contains(ons[i],ns))
						drawer.paintSubtree(offscreenGraphics, 
								getXScaleEq(), getYScaleEq(),
								destRect, ons[i], false, corr);

			//Repaint all subtrees that are newly selected
			if (ns != null)
				for (int i=0; i<ns.length; i++)
					//if (!contains(ns[i],ons))
						drawer.paintSubtree(offscreenGraphics, 
								getXScaleEq(), getYScaleEq(),
								destRect, ns[i], true, corr);
			*/

			drawer.paint(offscreenGraphics, 
					getXScaleEq(), getYScaleEq(),
					destRect, ns, corr);
		}
		repaint();
	}

	private void selectParent(TreeDrawerNode n)
	{
		TreeDrawerNode current = n;
		TreeDrawerNode parent = current.getParent();
		if (parent == null)
			return;
		/*
		if (current == parent.getLeft())
			current = parent.getRight();
		else
			current = parent.getLeft();
		drawer.paintSubtree(offscreenGraphics, 
				getXScaleEq(), getYScaleEq(),
				destRect, current, true);
		drawer.paintSingle(
				offscreenGraphics, getXScaleEq(), getYScaleEq(),
				destRect, current, true);
		redrawTree(new TreeDrawerNode[] {current}, new TreeDrawerNode[] {parent});
		repaint();
		*/
		geneSelection.setSelectedNode(parent.getId());
		geneSelection.notifyObservers();
	}

	private void selectRight(TreeDrawerNode n)
	{
		// XXX
		TreeDrawerNode right = n.getLeft();
		if (right.isLeaf()) return;
		/*
		drawer.paintSingle(offscreenGraphics,
				getXScaleEq(), getYScaleEq(), destRect, current, false);
		drawer.paintSubtree(offscreenGraphics, 
				getXScaleEq(), getYScaleEq(),
				destRect, current.getLeft(), false);
		redrawTree(new TreeDrawerNode[] {current}, new TreeDrawerNode[] {right});
		repaint();
				*/
		geneSelection.setSelectedNode(right.getId());
		geneSelection.notifyObservers();
	}

	private void selectLeft(TreeDrawerNode n)
	{
		// XXX
		TreeDrawerNode left = n.getRight();
		if (left.isLeaf()) return;
		/*
		drawer.paintSingle(offscreenGraphics,
				getXScaleEq(), getYScaleEq(), destRect, current, false);
		drawer.paintSubtree(offscreenGraphics, 
				getXScaleEq(), getYScaleEq(),
				destRect, current.getRight(), false);
		redrawTree(new TreeDrawerNode[] {current}, new TreeDrawerNode[] {left});
		repaint();
		*/
		geneSelection.setSelectedNode(left.getId());
		geneSelection.notifyObservers();
	}

	/**
	 * expect updates to come from map, geneSelection and drawer
	 */
	public void update(Observable o, Object arg) {
	   if(_debug) {
	      System.out.println(">>>>>> MultiSelectionGTRView: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
	   }
	   if (o == map) {
	      // System.out.println("Got an update from map");
	      offscreenValid = false;
	      repaint();
	   } else if (o == drawer) {
	      //System.out.println("Got an update from drawer");
	      offscreenValid = false;
	      repaint();
	   } else if (o == geneSelection) {
	      /*
	       * Not really bothered about this stuff right now
	       *
	       TreeDrawerNode cand = null;
	       if (geneSelection.getNSelectedIndexes() > 0) {
	      // This clause selects the array node if only a single array is selected.
	      if (geneSelection.getMinIndex() == geneSelection.getMaxIndex())
	      cand = drawer.getLeaf(geneSelection.getMinIndex());
	      // this clause selects the root node if all genes are selected.
	      else if ((geneSelection.getMinIndex() == map.getMinIndex()) 
	      && (geneSelection.getMaxIndex() == map.getMaxIndex()))
	      cand = drawer.getRootNode();
	      }
	      if ((cand != null)
	      && (cand.getId() != geneSelection.getSelectedNode()))
	      {
	      String id = cand.getId();
	      geneSelection.setSelectedNode(id);
	      geneSelection.notifyObservers();
	      }
	      else
	      {
	       */

	      /*
	       * New update process is as follows:
	       *
	       * 1. Convert list of node names from the gene selection
	       * to a list of TreeDrawerNodes
	       * 
	       * 2. Check that the selection of nodes has actually
	       * changed
	       *
	       * 3. If it has, update the selected indexes
	       *
	       * 4. Redraw the tree
	       *
	       * 5. Set old copy of node selection
	       *
	       * 6. Notify observers of the gene selection
	       */

	      // Stage 1.
	      String[] newHighlightedNames
		 = geneSelection.getHighlightedNodes().toArray(new String[0]);
	      TreeDrawerNode[] newHighlightedNodes = null;
	      if (newHighlightedNames != null)
	      {
		 newHighlightedNodes
		    = new TreeDrawerNode[newHighlightedNames.length];
		 for (int i=0; i<newHighlightedNames.length; i++)
		    newHighlightedNodes[i]
		       = drawer.getNodeById(newHighlightedNames[i]);
	      }

	      /*
		 String[] newSelectedNames
		 = geneSelection.getSelectedNodes().toArray(new String[0]);
		 TreeDrawerNode[] newSelectedNodes = null;
		 if (newSelectedNames != null)
		 {
		 newSelectedNodes
		 = new TreeDrawerNode[newSelectedNames.length];
		 for (int i=0; i<newSelectedNames.length; i++)
		 newSelectedNodes[i]
		 = drawer.getNodeById(newSelectedNames[i]);
		 }
	       */

	      /*if (MultiTreeSelection.hasReallyChanged(currentSelectedNodes, newSelectedNodes))
		{
	      // Stage 3.
	      //geneSelection.synchMap();
	      currentSelectedNodes = newSelectedNodes;
	      geneSelection.notifyObservers();
	      //printSelections(3);
	      }
	       */

	      redrawTree(newHighlightedNodes);

	      //printSelections(1);
	      // Stage 2.
	      /*
		 if (MultiTreeSelection.hasReallyChanged(currentHighlightedNodes, newHighlightedNodes))
		 {
	      // Stage 4.
	      //redrawTree(currentHighlightedNodes, newHighlightedNodes);
	      //printSelections(4);
	      // Stage 5.
	      currentHighlightedNodes = newHighlightedNodes;
	      //printSelections(5);
	      // Stage 6.
	      geneSelection.notifyObservers();
	      //printSelections(6);
	      }
	       */
	      //}
	   } else {
	      LogBuffer.println(viewName() + "Got an update from unknown " + o);
	   }
	   if(_debug) {
	      System.out.println("<<<<<< MultiSelectionGTRView: exit update, Observable = "+o.getClass().getSimpleName());
	   }
	} // update

	/*
	private void printSelections(int stage)
	{
		if (currentHighlightedNodes != null)
			for (int i=0; i<currentHighlightedNodes.length; i++)
				System.out.println("Old Node "+i+": "+currentHighlightedNodes[i].getId());
		if (geneSelection != null)
			geneSelection.printDetails(this.viewName()+" after stage "+stage);
	}
	*/

	private int hasChanged(TreeDrawerNode[] oldNodes, TreeDrawerNode[] newHighlightedNodes)
	{
		int changed = 0;
		if (oldNodes == null && newHighlightedNodes == null)
			changed = 0;
		else if (oldNodes == null)
			changed = 1;
		else if (newHighlightedNodes == null)
			changed = 2;
		else if (oldNodes.length != newHighlightedNodes.length)
			changed = 3;
		else
			for (int i=0; i<oldNodes.length; i++)
				if (!oldNodes[i].getId().equals(newHighlightedNodes[i].getId()))
					changed = 4;
		return changed;
	}

	// method from ModelView
	public String viewName() { return "MultiSelectionGTRView";}

	// method from ModelView
	public void updateBuffer(Graphics g) {
		//		System.out.println("GTRView updateBuffer() called offscreenChanged " + offscreenChanged + " valid " + offscreenValid + " yScaleEq " + getYScaleEq());
		if (offscreenChanged == true) offscreenValid = false;
		if ((offscreenValid == false) && (drawer != null)) {
			map.setAvailablePixels(offscreenSize.height);

			// clear the pallette...
			g.setColor(Color.white);
			g.fillRect
				(0,0, offscreenSize.width, offscreenSize.height);
			g.setColor(Color.black);

			//	calculate Scaling
			destRect.setBounds(0,0, offscreenSize.width,map.getUsedPixels());
			setXScaleEq( new LinearTransformation
					(drawer.getCorrMin(), destRect.x,
					 drawer.getCorrMax(), destRect.x + destRect.width));

			setYScaleEq(new LinearTransformation
					(map.getIndex(destRect.y), destRect.y,
					 map.getIndex(destRect.y + destRect.height), 
					 destRect.y + destRect.height));
			// System.out.println("yScaleEq " + getYScaleEq());
			// draw
			String[] names = geneSelection.getHighlightedNodes().toArray(new String[0]);
			TreeDrawerNode[] nodes = new TreeDrawerNode[names.length];
			for (int i=0; i<names.length; i++)
				nodes[i] = drawer.getNodeById(names[i]);
			drawer.paint(g, 
					getXScaleEq(), getYScaleEq(),
					destRect, nodes, geneSelection.getCorrelationValue());
		} else {
			//	        System.out.println("didn't update buffer: valid = " + offscreenValid + " drawer = " + drawer);
		}
	}


	// Mouse Listener 
	public void mouseClicked(MouseEvent e) {
	   if(_debug) {
	      System.out.println(">>>>>> MultiSelection GTRView: enter mouseClicked");
	   }
	   if (enclosingWindow().isActive() == false) {
	      if(_debug) {
		 System.out.println("mouseClicked returning, enclosing window not active");
	       }
	      return;
	   }
	   if(_debug) {
	      if (drawer == null) {
		 System.out.println("GTRView.mouseClicked() : drawer is null");
	      }
	      if (getXScaleEq() == null) {
		 System.out.println("GTRView.mouseClicked() : xscaleEq is null");
	      }
	   }
	   if ((drawer != null) && (getXScaleEq() != null)) {
	      // the trick is translating back to the normalized space...
	      TreeDrawerNode closest = drawer.getClosest (getYScaleEq().inverseTransform(e.getY()),
		       getXScaleEq().inverseTransform(e.getX()),
		       getXScaleEq().getSlope() / getYScaleEq().getSlope());

	      geneSelection.resetCorrelationValue();
	      int button = e.getButton();
	      if (button == MouseEvent.BUTTON1) {
		 if (e.isShiftDown()) {
		    if(_debug) {
		       System.out.println("calling geneSelection.setNode");
		    }
		    geneSelection.setNode(closest.getId(), !geneSelection.isSelected(closest.getId()));
		 } else {
		    if(_debug) {
		       System.out.println("calling geneSelection.setSelectedNode");
		    }
		    geneSelection.setSelectedNode(closest.getId());
		 }

		 //System.out.println("MultiSelection GTRView: mouseCLicked about to notifyObservers +++++++++++++++++++++");
		 geneSelection.notifyObservers((Object)"treeNodeClicked");

		 //Call a fake update so that the image stuff doesn't
		 //redraw before we've called synchMap().
		 //update((Observable)geneSelection, null);
	      }
	      /*
		 else if (SwingUtilities.isRightMouseButton(e))
		 {
		 geneSelection.setNode(closest.getId(), false);

		 geneSelection.resetCorrelationValue();
		 geneSelection.notifyObservers();
		 }
	       */
	   }
	   if(_debug) {
	      System.out.println("<<<<<< MultiSelection GTRView: exit mouseClicked");
	   }
	} // mouseClicked

	// method from KeyListener
	public void keyPressed(KeyEvent e)
	{
		String[] selectedNodes = geneSelection.getHighlightedNodes().toArray(new String[0]);
		TreeDrawerNode node = drawer.getNodeById(selectedNodes[0]);
		if (selectedNodes == null || selectedNodes.length != 1) {return;}
		int c = e.getKeyCode();	
		switch (c) {
			case KeyEvent.VK_UP:
				selectParent(node); break;
			case KeyEvent.VK_LEFT:
				if (node.isLeaf() == false)
					selectLeft(node);
				break;
			case KeyEvent.VK_RIGHT:
				if (node.isLeaf() == false)
					selectRight(node);
				break;
			case KeyEvent.VK_DOWN:
				if (node.isLeaf() == false)
				{
					// XXX
					TreeDrawerNode right = node.getLeft();
					TreeDrawerNode left = node.getRight();
					if (right.getRange() >  left.getRange())
						selectRight(node); 
					else
						selectLeft(node);
				}
				break;
		}
	}

	public void scrollToNode(String nodeName) {
		TreeDrawerNode node = drawer.getNodeById(nodeName);
		if (node != null) {
			int index = (int) node.getIndex();
			if (map.isVisible(index) == false) {
				map.scrollToIndex(index);
				map.notifyObservers();
			}
		}
	}
}
