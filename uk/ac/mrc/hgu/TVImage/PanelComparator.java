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
 * Used to sort image panels in a list.
 *
 * Currently, there's no way to sort a mixed collection of
 * nodes and genes according to their inorder position in the
 * tree.  Instead, nodes are placed first and genes second.
 *
 * TODO: investigate whether clever sorting of nodes and genes
 * is possible
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
class PanelComparator implements Comparator<ImageContainer>
{
	public PanelComparator() { super(); }

	public int compare(ImageContainer a, ImageContainer b)
	{
		int ai = a.getSortIndex();
		int bi = b.getSortIndex();
		if (ai < bi)
			return -1;
		else if (ai > bi)
			return 1;
		else
			return 0;
	}
}
