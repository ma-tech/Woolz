package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.util.*;
import javax.imageio.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.text.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/**
 * Displays information about the selected gene(s)
 *
 * @author Paul Smith <psmith@hgu.mrc.ac.uk>, Tom Perry
 * <tperry@hgu.mrc.ac.uk>
 *
 * TODO: Currently won't display filename or date of database
 * extraction.  Also won't resize consistently, small number of
 * genes = large buttons. large number of genes = small buttons
 */
class SummaryView extends JScrollPane implements Observer, ComponentListener {

	private ImageViewManager tvim;
	//private HeaderInfo gtrInfo = null;
	private Collection<String> selectedNodes = new HashSet<String>();
	private Collection<String> selectedGenes = new HashSet<String>();

	private Set<String> requests = new HashSet<String>();

	private String description = null;
	private XmlAppSettings appSettings;

	//private JScrollPane scrollPane;
	private JEditorPane textComponent = new JEditorPane("text/html","");

	private boolean _debug = false;

	public SummaryView(ImageViewManager m)
	{
		super();
		tvim = m;

		//scrollPane = new JScrollPane();
		//add(scrollPane);
		description = "Summary";
		addComponentListener(this);
		textComponent.setEditable(false);
		//scrollPane.setViewportView(textComponent);
		setViewportView(textComponent);
		textComponent.addHyperlinkListener(new HyperlinkListener()
		{
			public void hyperlinkUpdate(HyperlinkEvent e)
			{
				if (e.getEventType().equals(HyperlinkEvent.EventType.ACTIVATED))
					tvim.displayURL(e.getURL());
			}
		});
		
		setPreferredSize(new Dimension(getVisibleRect().width, getVisibleRect().height));

		appSettings = tvim.getAppSettingsObject();
	}

	public void componentResized(ComponentEvent e) { doMediaLayout(); }
	public void componentHidden(ComponentEvent e) {}
	public void componentMoved(ComponentEvent e) {}
	public void componentShown(ComponentEvent e) {}

	public String viewName() { return "ImageView";}
	public String viewDescription() { return description; }

	public ImageViewManager getManager() { return tvim; }

	/*
	private JPanel createNodePanel()
	{
		JPanel jp = new JPanel();

		jp.setLayout(new GridLayout(5, 1));

		if ((selectedNodes.isEmpty()) && !(selectedGenes.isEmpty()))
		{
			jp.add(createRequiredJLabel("No node is selected", false));
		}
		else if (selectedNodes.size() > 1) 
		{
			jp.add(createRequiredJLabel("Multiple nodes are selected", false));
		}
		else
		{
			for (String s : selectedNodes)
			{
				jp.add(createRequiredJLabel("NODE ID: " + s, false));
				jp.add(createRequiredJLabel("CORRELATION: " + tvim.getCorrelation(s), false));
				jp.add(createRequiredJLabel("NUMBER OF LEAVES: " + selectedGenes.size(), false));
				jp.add(createRequiredJLabel("DATASET: " + "", false));
				jp.add(createRequiredJLabel("FILENAME: ", false));
			}
		}

		jp.setVisible(true);
		jp.setPreferredSize(new Dimension(getVisibleRect().width,
					getPreferredSize().height));
		return jp;
	}
	*/

	private void setupTextPanel()
	{
		String text = "";
		text += "<style type=\"text/css\">";
		text += "table { border-collapse: collapse; }";
		//text += "font: 120% arial, helvetica, sans-serif;";
		text += "</style>";

		String[] attributes = new String[] {"Node ID","Correlation",
			"Number of leaves","Dataset","Date of clustering"};

		String[] values = new String[attributes.length];
		if (selectedNodes.size() == 1)
		{
			String node = selectedNodes.iterator().next();
			values[0] = node;
			values[1] = tvim.getCorrelation(node);
		}
		values[2] = String.valueOf(selectedGenes.size());

		text += "<table border=\"1\">";
		for (int i=0; i<attributes.length; i++)
		{
			text += "<tr><th align=\"left\">";
			text += attributes[i];
			text += "</th><td>";
			if (values[i] != null)
				text += values[i];
			text += "</td></tr>";
		}
		text += "</table>";

		if (!selectedGenes.isEmpty())
		{
			boolean showGenes = true;
			boolean showEMAGE = true;

			text += "<table border=\"1\"><tr>";
			if (showGenes)
				text += "<th align=\"left\">Genes</th>";
			if (showEMAGE)
				text += "<th align=\"left\">EMAGE ID</th>";
			text += "</tr>";
			for (String s : selectedGenes)
			{
				String geneName = tvim.getGeneName(s);
				String emageId = tvim.getEMAGEID(s);
				text += "<tr>";
				if (showGenes)
				{
					text += "<td><a href=\"";
					text += appSettings.getGeneNameQueryString(geneName);
					text += "\">";
					text += geneName;
					text += "</a></td>";
				}
				if (showEMAGE)
				{
					String emageNum;
					if (emageId.length() > 6) 
						emageNum = emageId.substring(6);
					else
						emageNum = emageId;
					text += "<td><a href=\"";
					text += appSettings.getEmagePageLocation(emageNum);
					text += "\">";
					text += emageId;
					text += "</a></td>";
				}
				text += "</tr>";
			}
			text += "</table>";
		}
		textComponent.setText(text);
		//textComponent.setVisible(true);
		//textComponent.setPreferredSize(new Dimension(getVisibleRect().width,
					//getPreferredSize().height));
	}

	private synchronized void doMediaLayout()
	{
		//removeAll();
		//scrollPane.setPreferredSize(new Dimension(
					//getVisibleRect().width,
					//getVisibleRect().height));
		setPreferredSize(new Dimension(
					getVisibleRect().width,
					getVisibleRect().height));
		setupTextPanel();
		revalidate();
		repaint();
	}

	/** 
	 * Method from java.util.Observer that is called when an
	 * observed object is changed in some way.
	 *
	 * The parent ImageViewManager is an observer of its
	 * geneSelection object, and this class becomes involved
	 * because any calls to update() in the ImageViewManager class
	 * are passed to all of the Panels that it manages.
	 *
	 * @param o 
	 * @param arg 
	 */
	public synchronized void update(Observable o, Object arg) {

	   if(_debug) {
	      System.out.println(">>>>>> SummaryView: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
	   }
	   if (o.equals(tvim.getGeneSelection())) {
	      Collection<String> newSelectedNodes = tvim.getGeneSelection().getSelectedNodes();
	      Collection<String> newSelectedGenes = tvim.getSelectedGenes();

	      //want to setup details for the node selected
	      selectedNodes = newSelectedNodes; 
	      selectedGenes = newSelectedGenes;
	   }
	   /*
	    * Assume that any other updates must be from the image
	    * server, though this may turn out to be a bad assumption.
	    *
	    * Some checking may be a good idea, but this requires
	    * identification of ImageLoaders which are currently
	    * anonymous (see ImageServer.java).
	    */
	   doMediaLayout();
	   if(_debug) {
	      System.out.println("<<<<<< SummaryView: exit update, Observable = "+o.getClass().getSimpleName());
	   }
	}
}
