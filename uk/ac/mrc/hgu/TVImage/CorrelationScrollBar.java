/* file name  : CorrelationScrollBar.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
 * created    : Mon 03 Sep 2007 09:34:33 BST
 */
package uk.ac.mrc.hgu.TVImage;

import edu.stanford.genetics.treeview.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.util.*;

/** 
 * Allows the user to select a correlation value by dragging a
 * slider.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
public class CorrelationScrollBar extends JSlider implements Observer {

	private ImageViewManager tvim;
	private static final int MAX = 100;
	private static final boolean RESPONSIVE_SCROLLING = true;
	private MultiTreeSelectionI geneSelection;
	private boolean _debug;

	public CorrelationScrollBar(ImageViewManager m) {

	   super(JSlider.HORIZONTAL,0,MAX,MAX/2);
	   tvim = m;

	   addMouseListener(new MouseListener() {
	      public void mouseReleased(MouseEvent me) {
		 if (!RESPONSIVE_SCROLLING) {
		    double corr = ((double)getValue())/((double)MAX);
		    tvim.selectNodesByCorrelation(corr);
		 }
	      }
	      public void mouseClicked(MouseEvent me) {}
	      public void mousePressed(MouseEvent me) {}
	      public void mouseEntered(MouseEvent me) {}
	      public void mouseExited(MouseEvent me) {}
	   });

	   addChangeListener(new ChangeListener() {
	      public void stateChanged(ChangeEvent e) {
		 double corr = ((double)getValue())/((double)MAX);

		 //tvim.setCorrelationView(true);

		 //System.out.println("CorrelationScrollBar: stateChanged RESPONSIVE_SCROLLING = "+RESPONSIVE_SCROLLING);
		 if (RESPONSIVE_SCROLLING) {
		    tvim.selectNodesByCorrelation(corr);
		 } else {
		    tvim.getGeneSelection().setCorrelationValue(corr);
		    tvim.getGeneSelection().deleteObserver(tvim);
		    tvim.getGeneSelection().notifyObservers();
		    tvim.getGeneSelection().addObserver(tvim);
		 }
	      }
	   });
	}

	public void update(Observable o, Object arg) {
	   if(_debug) {
	      System.out.println(">>>>>> CorrelationScrollBar: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
	   }
	   if (o == tvim.getGeneSelection()) {
	      int value = (int)(geneSelection.getCorrelationValue()*MAX);
	      if (value >= 0 && value <= MAX) {
		 setValue(value);
	      }
	   }
	   if(_debug) {
	      System.out.println("<<<<<< CorrelationScrollBar: exit update, Observable = "+o.getClass().getSimpleName());
	   }
	}
	
	/* 
	 * Putting this method here is not ideal but this class needs
	 * to observe some geneSelection object and this appears to be
	 * the only way to do it.
	 */
	public void setGeneSelection(TreeSelectionI geneSelection) {
	   if (this.geneSelection != null) {
	      this.geneSelection.deleteObserver(this);
	   }
	   this.geneSelection = (MultiTreeSelectionI)geneSelection;
	   update((Observable)geneSelection, null);
	   if (this.geneSelection != null) {
	      this.geneSelection.addObserver(this);
	   }
	}
}
