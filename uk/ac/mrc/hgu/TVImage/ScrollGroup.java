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

class ScrollGroup extends HashMap<Integer,JScrollBar>
	implements AdjustmentListener, MouseListener,
						 MouseWheelListener
{
	ScrollGroup() { super(); }
	JScrollBar sourceScroll = null;

	public void adjustmentValueChanged(AdjustmentEvent e)
	{
		if (sourceScroll == null)
			return;

		for (JScrollBar jsb : this.values())
			if (!jsb.equals(sourceScroll))
			{
				int max = sourceScroll.getMaximum();
				if (max == 0)
					return;
				else
				{
					jsb.setValue(sourceScroll.getValue()
							* jsb.getMaximum() / max);
				}
			}
	}

	public void mousePressed(MouseEvent e)
	{ sourceScroll = (JScrollBar)e.getSource(); }

	public void mouseReleased(MouseEvent e)
	{ sourceScroll = null; }

	public void mouseClicked(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}

	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int count = e.getWheelRotation();
		Object src = e.getSource();

		if (src instanceof JScrollBar)
			sourceScroll = (JScrollBar)src;
		else if (src instanceof JScrollPane)
			sourceScroll = ((JScrollPane)src).getVerticalScrollBar();

		sourceScroll.setValue(sourceScroll.getValue() + count*sourceScroll.getMaximum()/20);
	}
}
