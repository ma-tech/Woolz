package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.util.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

public class MultiSelectionLeftTreeDrawer extends LeftTreeDrawer {
	// Holds previous correlation value so that old line may be
	// painted over
	private double currentCorr = -1;

	public void paint(Graphics graphics, LinearTransformation xScaleEq, 
			LinearTransformation yScaleEq, Rectangle dest, 
			TreeDrawerNode[] selected, double corr) {	
		if ((getRootNode() == null) || (getRootNode().isLeaf() == true))
			System.out.println("Root node is null or leaf!");
		else {
			// recursively drawtree...
			NodeDrawer nd = new NodeDrawer(
					graphics, xScaleEq, yScaleEq, 
					selected, dest, corr);
			nd.draw(getRootNode());
		}
	}

	public void paintSubtree(Graphics graphics, LinearTransformation xScaleEq, 
			LinearTransformation yScaleEq, Rectangle dest, 
			TreeDrawerNode root, boolean isSelected, double corr) {	
		if ((root == null) || (root.isLeaf() == true))
			return;
		else {
			if (yScaleEq == null) {
				LogBuffer.println("yScaleEq was null in LeftTreeDrawer.paintSubTree!");
				Exception e = new Exception();
				e.printStackTrace();
			}
			// recursively drawtree...
			NodeDrawer nd = new NodeDrawer(
					graphics, xScaleEq, yScaleEq, 
					null, dest, corr);
			nd.isSelected = isSelected;
			nd.draw(root);
		}
	}

	public void paintSubtree(Graphics graphics, LinearTransformation xScaleEq, 
			LinearTransformation yScaleEq, Rectangle dest, 
			TreeDrawerNode root, 
			TreeDrawerNode[] selected, double corr) {	
		if ((root == null) || (root.isLeaf() == true))
			return;
		else {
			// recursively drawtree...
			NodeDrawer nd = new NodeDrawer(
					graphics, xScaleEq, yScaleEq, 
					selected, dest, corr);
			nd.draw(root);
		}
	}

	public void paintSingle(Graphics graphics, LinearTransformation xScaleEq, 
			LinearTransformation yScaleEq, Rectangle dest, 
			TreeDrawerNode root, boolean isSelected, double corr) {	
		if ((root == null) || (root.isLeaf() == true))
			return;
		else {
			// just draw single..
			NodeDrawer nd = new NodeDrawer(
					graphics, xScaleEq, yScaleEq, 
					null, dest, corr);
			nd.isSelected = isSelected;
			if (root.isLeaf() == false)
				nd.drawSingle(root);
			else
				System.err.println("Root was leaf?");
		}
	}


	/**
	 * this is an internal helper class which does a sort of recursive drawing
	 * that's actually implemented with iteration.
	 * 
	 * @author Alok Saldanha <alok@genome.stanford.edu>
	 * @version Alpha
	 */
	class NodeDrawer
	{
		boolean isSelected = false;
		private Color sel_color = Color.red;
		private Color node_color = Color.blue;
		private Graphics graphics;
		private TreeDrawerNode[] selected;
		private LinearTransformation xT, yT;

		private double minInd, maxInd;
		private Rectangle dest;
		private double corr;

		/**
		 * The constructor sets the variables
		 * 
		 * @param g         The graphics object to print to
		 * @param xScaleEq The equation to be applied to scale the index of the nodes to graphics object
		 * @param yScaleEq The equation to be applied to scale the correlation of the nodes to the graphics object
		 * maybe foreground color, selection color and node color should be options?
		 */
		public NodeDrawer(Graphics g, LinearTransformation xScaleEq, 
				LinearTransformation yScaleEq, TreeDrawerNode[] sel,
				Rectangle d, double c) 
		{
			if (yScaleEq == null) {
				LogBuffer.println("yScaleEq was null!");
				return;
			}
			graphics = g;
			selected = sel;
			xT = xScaleEq;
			yT = yScaleEq;
			dest = d;
			corr = c;
			if (dest != null) {
				minInd = (int) yScaleEq.inverseTransform(dest.y);
				maxInd = (int) yScaleEq.inverseTransform(dest.y + dest.height) + 
					1;
			}
		}

		private boolean contains(TreeDrawerNode n, TreeDrawerNode[] ns)
		{
			if (n == null || ns == null)
				return false;
			for (int i=0; i<ns.length; i++)
				if (n == ns[i])
					return true;
			return false;
		}

		/** 
		 * the draw method actually does the drawing
		 */
		public void draw(TreeDrawerNode startNode)
		{
			graphics.setColor(Color.white);
			graphics.drawLine((int)(dest.width*currentCorr), 0, (int)(dest.width*currentCorr), dest.height);
			graphics.setColor(Color.red);
			graphics.drawLine((int)(dest.width*corr), 0, (int)(dest.width*corr), dest.height);
			graphics.setColor(Color.black);
			Stack<TreeDrawerNode> remaining = new Stack<TreeDrawerNode>();
			remaining.push(startNode);
			while (!remaining.empty())
			{
				TreeDrawerNode node = (TreeDrawerNode) remaining.pop();
				// just return if no subkids visible.
				if ((node.getMaxIndex() < minInd) ||
						(node.getMinIndex() > maxInd))
					continue;

				// handle selection...
				if (contains(node,selected))
				{
					if (isSelected == false)
					{
						isSelected = true;
						// push onto stack, so we know when we're finished with the
						// selected subtree..
						remaining.push(node);
					}
					else
					{
						// isSelected is true, so we're pulling the selected node off
						// the second time.
						//
						// additional note from Tom:
						// children of selected nodes can now be selected as well, so we
						// need to make sure that any selected descendants of selected
						// nodes don't unset the isSelected-ness here
						//
						// TODO: this implementation is inefficient - fix?

						TreeDrawerNode currentDescendant = node;
						boolean selectedParent = false;
						while (currentDescendant != null)
						{
							currentDescendant = currentDescendant.getParent();
							if (contains(currentDescendant, selected))
								selectedParent = true;
						}
						if (!selectedParent)
						{
							isSelected = false;
							continue;
						}
					}
				}
				// lots of stack allocation...
				TreeDrawerNode left = node.getLeft();
				TreeDrawerNode right = node.getRight();
				if (!left.isLeaf())
					remaining.push(left);
				if (!right.isLeaf())
					remaining.push(right);
				// finally draw
				drawSingle(node);
			}

			if (corr != currentCorr)
				currentCorr = corr;
		}

		private void drawSingle(TreeDrawerNode node) {
			TreeDrawerNode left = node.getLeft();
			TreeDrawerNode right = node.getRight();
			if (xT == null) 
				System.err.println("xt was null");
			if (right == null)
				System.err.println("right was null");
			int rx = (int) xT.transform(right.getCorr());
			int lx = (int) xT.transform(left.getCorr());
			int tx = (int) xT.transform(node.getCorr());

			int ry = (int) yT.transform(right.getIndex() + .5);
			int ly = (int) yT.transform(left.getIndex() + .5);
			int ty = (int) yT.transform(node.getIndex() + .5);

			// System.out.println("rx = " + rx + ", ry = " + ry + ", lx = " + lx + ", ly = " + ly + "corr " + node.getCorr());

			// oval first?...
			//	    graphics.setColor(node_color);
			//	    graphics.drawOval(tx - 1,ty - 1,2,2);

			//draw our (flipped) polyline...
			Random r = new Random();
			if (isSelected) 
				//graphics.setColor(new Color(r.nextInt(127),r.nextInt(127),r.nextInt(127)));
				graphics.setColor(sel_color);
			else
				//graphics.setColor(new Color(r.nextInt(127)+128,r.nextInt(127)+128,r.nextInt(127)+128));
				graphics.setColor(node.getColor());

			//graphics.drawPolyline(new int[] {rx, tx, tx, lx},
			//		new int[] {ry, ry, ly, ly}, 4);

			graphics.drawLine(tx,ry,tx,ly);

			if (corr >= 0 && right.isLeaf())
				graphics.setColor(sel_color);
			graphics.drawLine(rx,ry,tx,ry);
			if (corr >= 0 && left.isLeaf())
				graphics.setColor(sel_color);
			graphics.drawLine(tx,ly,lx,ly);
		}
	}
}
