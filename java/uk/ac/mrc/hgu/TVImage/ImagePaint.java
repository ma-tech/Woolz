package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.util.*;
import javax.imageio.*;
import java.io.*;
import java.text.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/**
 * TODO: write description here
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 * @version 
 */
class ImagePaint extends ModelView implements MouseMotionListener {

	private BufferedImage canvas;
	private BufferedImage model;
	private TVWlzObject comparisonImg;
	private Graphics2D canvasGraphics;
	private Dimension modelDimensions;
	private Point mousePosition = new Point();
	private ImageViewManager tvim;
	private ImageServer imgProducer;
	private boolean _debug = false;

	//TODO: We need some way of generating this on demand
	private static final String modelPath = "/net/homehost/export/home/tperry/models/";

	/** 
	 * 
	 * 
	 * @param tvim 
	 */
	public ImagePaint(ImageViewManager tvim)
	{
		super();
		panel = this;
		this.tvim = tvim;
		addMouseListener(this);
		addMouseMotionListener(this);

		try
		{
			model = ImageIO.read(new File(modelPath + "model.jpg"));
		}
		catch (IOException ioe) { ioe.printStackTrace(); }

		modelDimensions = new Dimension(model.getWidth(null),model.getHeight(null));
    
		canvas = new BufferedImage(modelDimensions.width, modelDimensions.height, BufferedImage.TYPE_BYTE_BINARY);
		canvasGraphics = (Graphics2D)canvas.createGraphics();

		//imgProducer.loadImagesMT(Collections.singleton("GENE9X"),(Observer)this,ImageViewManager.TYPE_HEATMAP);
	}

	public String viewName() { return "ImagePaint";}
	public String viewDescription() { return "Painter";}

	public void mousePressed(MouseEvent me)
	{
		mousePosition.x = me.getX();
		mousePosition.y = me.getY();
		canvasGraphics.setColor(Color.white);
		int d = 10;
		canvasGraphics.fillOval(mousePosition.x-d/2,mousePosition.y-d/2,d,d);
		repaint();
	}

	public void mouseDragged(MouseEvent me)
	{
		int nx = me.getX();
		int ny = me.getY();
		canvasGraphics.setColor(Color.white);
		canvasGraphics.setStroke(new BasicStroke((float)10.0,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND));
		canvasGraphics.drawLine(mousePosition.x,mousePosition.y,nx,ny);
		mousePosition.x = nx;
		mousePosition.y = ny;
		repaint();
	}

	public void updateBuffer(Graphics g)
	{
		Graphics2D g2 = (Graphics2D)g;
		g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
		g2.drawImage(model, 0, 0, this);
		g2.drawImage(canvas, 0, 0, this);

		g2.setColor(Color.black);

		double correlation = TVWlzObject.getJaccardCorrelation(comparisonImg,new TVWlzObject(canvas));
		g2.drawString("Correlation: " + correlation,0,400);
	}

	/**
	 * Callback for image loader
	 *
	 * @param o   the java.util.Observable object under observation
	 * @param arg optional argument
	 */
	public void update(Observable o, Object arg) {
	   if(_debug) {
	      System.out.println(">>>>>> ImagePaint: enter update, Observable = "+o.getClass().getSimpleName()+", arg = "+arg);
	   }
	   //if ((String)arg == "GENE9X")
	   //	comparisonImg = new TVWlzObject(tvim.getMedia(Collections.singleton("GENE9X"),ImageViewManager.TYPE_HEATMAP).get(0).getImage());
	   //else
	   repaint();
	   if(_debug) {
	      System.out.println("<<<<<< ImagePaint: exit update, Observable = "+o.getClass().getSimpleName());
	   }
	}
}
