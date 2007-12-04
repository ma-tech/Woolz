/* file name  : CorrelationBox.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
 * created    : Mon 03 Sep 2007 09:17:07 BST
 */
package uk.ac.mrc.hgu.TVImage;

import edu.stanford.genetics.treeview.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

/** 
 * Shows the current correlation value held by the
 * MultiTreeSelection object.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
public class CorrelationBox extends JPanel implements Observer, TVTypes {

	private ImageViewManager tvim;
	private MultiTreeSelectionI geneSelection;
	private JTextField box;
	private JLabel label;
	private boolean _debug = false;

	public CorrelationBox(ImageViewManager m, int height) {
	   super();
	   tvim = m;
	   box = new JTextField(4);
	   String text = "Correlation: ";
	   label = new JLabel(text);
	   add(label);
	   add(box);

	   FlowLayout lm = new FlowLayout(FlowLayout.LEADING);
	   lm.setHgap(0);
	   lm.setVgap(0);

	   setBackground(TVTypes.BGCOLOR);
	   setPreferredSize(new Dimension(getVisibleRect().width, height));
	   setLayout(lm);

	   box.addActionListener(new ActionListener() {
	      public void actionPerformed(ActionEvent ae) {
		 //tvim.setCorrelationView(true);
		 tvim.selectNodesByCorrelation(Double.parseDouble(box.getText()));
	      }
	   });
	}

	public void update(Observable o, Object arg) {
	   if(_debug) {
	      System.out.println(">>>>>> CorrelationBox: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
	   }
	   double value = geneSelection.getCorrelationValue();
	   if (value >= 0 && value <= 1) {
	      box.setText(String.valueOf(value));
	   } else {
	      box.setText("");
	   }
	   if(_debug) {
	      System.out.println("<<<<<< CorrelationBox: exit update, Observable = "+o.getClass().getSimpleName());
	   }
	}

	public void setGeneSelection(TreeSelectionI geneSelection)
	{
		if (this.geneSelection != null) this.geneSelection.deleteObserver(this);
		this.geneSelection = (MultiTreeSelectionI)geneSelection;
		update((Observable)geneSelection, null);
		if (this.geneSelection != null) this.geneSelection.addObserver(this);
	}
}
