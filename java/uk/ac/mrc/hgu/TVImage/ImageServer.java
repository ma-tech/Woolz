/* file name  : ImageServer.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>, Paul Smith <psmith@hgu.mrc.ac.uk>
 * created    : Fri 13 Jul 2007 15:12:34 BST
 * copyright  : 
 *
 * modifications:
 * Paul Smith - Modified to load RAW and ANNOTATED genes via http (and ignore nodes for these types of images)
 */

package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.imageio.*;
import java.lang.reflect.Method;
import java.lang.reflect.InvocationTargetException;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

/** 
 * Handles everything to do with the retrieval of images.
 *
 * Currently, image loads are multithreaded, while image getting
 * is single-threaded.  This makes sense because we want to limit
 * mutithreaded activity to this class for the time being.
 *
 * @author Tom Perry, Paul Smith
 */
abstract class ImageServer {

	private Map<String,BufferedImage> imageMap;
	private static long time = 0;

	protected ImageViewManager tvim = null;
	protected XmlAppSettings xml = null;

	protected static final boolean IMAGE_CACHING = false;

	private boolean _debug = false;

//---------------------------------------------------------------------------
	ImageServer(ImageViewManager tvim, XmlAppSettings xml) {
	   this.tvim = tvim;
	   this.xml = xml;
	   // Set imageMap capacity to known tree size?
	   //imageMap = new HashMap<String,BufferedImage>(
	   //tvim.getNumNodes()*2, (float)1.0);
	   // Or let java work it out?
	   imageMap = new HashMap<String,BufferedImage>();
	}
	
//---------------------------------------------------------------------------
	/** 
	 * Gets an image from the hash table.  Single-threaded.
	 *
	 * Images may only be requested in a packaged form provided by
	 * getMedia() so this method is protected.
	 * 
	 * @param id the requested node ID
	 * @return the BufferedImage corresponding to the id, or null
	 * if no such image exists in the hash table
	 */
	protected BufferedImage getImage(String id) {
	   BufferedImage img = (BufferedImage)imageMap.get(id);
	   assert img != null : "Failed to get image for " + id;
	   if(img == null) {
	      //System.out.println("FAILED TO GET "+getClass().getSimpleName()+": "+id+" IMAGE FROM IMAGE MAP ((((((((((((((((((((((((((((");
	   } else {
	      //System.out.println("GOT "+getClass().getSimpleName()+": "+id+" IMAGE FROM IMAGE MAP ))))))))))))))))))))))))))");
	   }
	   return img;
	}

//---------------------------------------------------------------------------
	/** 
	 * Provides the path to the image for the requested node.
	 * 
	 * @param nodeID requested node
	 * @return path to requested node
	 */
	protected abstract String getImagePath(String nodeID);

//---------------------------------------------------------------------------
	/** 
	 * Gets the ImageContainer for the requested node or gene ID.
	 * 
	 * @param Set<String>mediaIDs requested node or gene ID
	 * @return object extending ImageContainer
	 */
	public abstract ImageContainer getMedia(String id);

//---------------------------------------------------------------------------
	/** 
	 * Checks whether an image that corresponds to the provided
	 * node ID exists in the hash table.
	 * 
	 * @param id the requested image identifier
	 * @return true if the hash table returns a non-null image
	 */
	protected boolean isLoaded(String id) {
	   //System.out.println(">>>>>> ImageServer: enter / exit isLoaded("+id+")");
	   return (imageMap.get(id) != null);
	}

//---------------------------------------------------------------------------
	/** 
	 * Load an image multi-threadedly.
	 *
	 * @param id the requested node or gene ID, in the form
	 * {NODE,GENE}??X
	 *
	 * @param requester a reference to the class instance that
	 * requested the image.  A reference is preferable because
	 * there may be numerous instances of a particular class that
	 * are all able to load images.
	 */
	public void loadImagesMT(Set<String> nodeIDs, Observer requester) {
	   for (String s : nodeIDs) {
	      new Thread(new ImageLoader(s, requester)).start();
	   }
	}

//---------------------------------------------------------------------------
	public void loadImagesST(Set<String> nodeIDs, ImageView requester, boolean isTop) {
	   //_debug = true;
	   if(_debug) {
	      System.out.println(">>>>>> ImageServer: enter loadImagesST, requester = "+requester.viewDescription()+" isTop = "+isTop);
	      for(String id : nodeIDs) {
		 System.out.println(id);
	      }
	   }
	   for (String nodeID : nodeIDs) {
	      if (!isLoaded(nodeID)) {
		 //System.out.println(nodeID +" not yet loaded");
		 startClock();
		 BufferedImage img = fetchImage(nodeID);
		 ImageViewManager.printDebugMessage(
		       "Loaded an image for "+nodeID+" in "+stopClock()+"ms");
		 assert img != null;
		 if(_debug) {
		    if(img == null) {
		       System.out.println("PUT NULL "+requester.viewDescription()+" IMAGE "+nodeID+" IN IMAGE MAP (((((((((((");
		    } else {
		       System.out.println("SUCCESSFULLY PUT "+requester.viewDescription()+" IMAGE "+nodeID+" IN IMAGE MAP ))))))))");
		    }
		 }
		 BufferedImage scaledImage = scaleImage(img);
		 imageMap.put(nodeID, scaledImage);
		 assert isLoaded(nodeID);
		 img = null;
	      } else {
		 if(_debug) {
		    System.out.println("!!!!!!!!!!!!!!!!!!!!!!  nodeID "+nodeID+" already loaded !!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		 }
		 ImageViewManager.printDebugMessage("Loading an image for "+nodeID+"... already loaded.");
	      }
	      //System.out.println("updating requester with "+nodeID);
	      //requester.update(null, nodeID);
	   }

	   requester.addLoadedImages(nodeIDs, this, isTop);

	   if(_debug) {
	      System.out.println("<<<<<< ImageServer: exit loadImagesST, requester = "+requester.viewDescription());
	   }
	   //_debug = false;
	}

//---------------------------------------------------------------------------
	/** 
	 * Wait for a random period.
	 */
	private synchronized void waitForAWhile() {
	   waitForAWhile(100+(new java.util.Random()).nextInt(500));
	}

//---------------------------------------------------------------------------
	/** 
	 * Attempts to simulate network latency by forcing a wait for
	 * the specified number of milliseconds.  Only one thread may
	 * be 'waiting' at a time.  Other threads will have to wait for
	 * an opportunity to 'wait', if that makes sense...
	 * 
	 * @param ms number of milliseconds to wait for
	 */
	private synchronized void waitForAWhile(int ms) {
	   synchronized(imageMap) {
	      try {
		 imageMap.wait(ms);
	      }
	      catch (InterruptedException ie) {}
	   }
	}

//---------------------------------------------------------------------------
	/** 
	 * Gets an image for the specified node and stores it in a map.
	 * Images may be loaded from a local file or from the web.
	 * 
	 * @param nodeID node for which an image should be fetched
	 */
	private BufferedImage fetchImage(String nodeID) {
	   //System.out.println(">>>>>> ImageServer: enter fetchImage "+nodeID);
	   String path = "";
	   BufferedImage img = null;

	   // Callers should check whether the image is loaded before
	   // asking for it to be fetched from the web.
	   assert !isLoaded(nodeID);

	   path = getImagePath(nodeID);
	   assert path != null && path != "";
	   //System.out.println("ImageServer.fetchImage: path = "+path);

	   if (path.startsWith("http")) {
	      try {
		 java.net.URL imgURL = new java.net.URL(path); 
		 img = ImageIO.read(imgURL); 
	      }
	      catch (MalformedURLException mue) {
		 System.err.println("Malformed URL: " + path);
	      }
	      catch (IOException ioe) {
		 System.err.println("Failed to load image from " + path);
	      }
	   } else {
	      try {
		 img = ImageIO.read(new File(path));
	      }
	      catch (IOException ioe) {
		 System.err.println("Failed to load image from " + path); }
	   }
	   //waitForAWhile(200);
	   assert img != null;
	   if(img == null) {
	      //System.out.println("<<<<<< ImageServer: exit fetchImage "+nodeID+" FAIL");
	   } else {
	      //System.out.println("<<<<<< ImageServer: exit fetchImage "+nodeID+" SUCCESS");
	   }
	   return img;
	} // fetchImage

//---------------------------------------------------------------------------
	/**
	 *   Removes BufferedImage from imageMap.
	 *   @param id the key for the image to be removed.
	 *   @param removeKey if true remove the key-value pair from the HashMap (otherwise set the value to null).
	 *   @return true if the image was non-null.
	 */
	protected boolean removeImage(String id, boolean removeKey) {

	   boolean ret = false;
	   BufferedImage img = null;

	   if(imageMap.containsKey(id)) {
	      img = imageMap.get(id);
	      if(img != null) {
		 ret = true;
		 img = null;
	      }
	      imageMap.remove(id);
	      if(!removeKey) {
		 imageMap.put(id, null);
	      }
	   }
	   return ret;
	}

//---------------------------------------------------------------------------
   /**
    *   Remove all entries from imageMap.
    */
    protected void clearAll() {

       BufferedImage img = null;
       Set<String> keys = imageMap.keySet();
       for(String str : keys) {
          img = imageMap.get(str);
	  if(img != null) {
	     img = null;
	  }
       }
       imageMap.clear();
    }
//---------------------------------------------------------------------------
   /**
    *   Makes a scaled BufferedImage.
    *   (based on comp.lang.java.gui reply by Shannon Hickey - Swing Team)
    *   Google:	scaling a BufferedImage w/ Java 2D
    *   @param bi the unscaled image.
    *   @return the scaled image.
    */
   private BufferedImage scaleImage(BufferedImage bi) {

      int H = bi.getHeight();
      int W = bi.getWidth();
      int scaledImageH = ImageContainer.imagePanelH;
      int scaledImageW = ImageContainer.imagePanelW;

      if(H <= scaledImageH && W <= scaledImageW) {
	 // no need to scale image to fit container
	 //System.out.println("ImageContainer.scaleImg: returning, img fits container already");
         return bi;
      }

      double ratio = ((double)H/(double)W);
      boolean portrait = (ratio >= 1.0); // the image is taller than it is wide (or square).
      //System.out.println("Image H/W = "+ratio);

      int scaledW = 0;
      int scaledH = 0;
      double newRatio = 1.0;

      if(portrait) {
	 //System.out.println("portrait");
	 scaledH = scaledImageH;
	 scaledW = (int)(Math.ceil((scaledImageH / ratio)));
	 //System.out.println("scaled img is "+scaledW+" wide X "+scaledH+" high");
	 newRatio = ((double)scaledW/(double)scaledImageW);
	 //System.out.println("new ratio = "+ratio);
	 if(newRatio > 1.0) {
	    scaledH /= newRatio;
	    scaledW /= newRatio;
	 }
	 //System.out.println("scaled img is now "+scaledW+" wide X "+scaledH+" high");
      } else {
	 //System.out.println("landscape");
	 scaledW = scaledImageW;
	 scaledH = (int)(Math.ceil((scaledImageW * ratio)));
	 //System.out.println("scaled img is "+scaledW+" wide X "+scaledH+" high");
	 newRatio = ((double)scaledH/(double)scaledImageH);
	 //System.out.println("new ratio = "+ratio);
	 if(newRatio > 1.0) {
	    scaledH /= newRatio;
	    scaledW /= newRatio;
	 }
	 //System.out.println("scaled img is now "+scaledW+" wide X "+scaledH+" high");
      }

      BufferedImage scaledImg = new BufferedImage(scaledW, scaledH, bi.getType());
      Graphics2D g = (Graphics2D)(scaledImg.getGraphics());

      // Frame implements imageObserver interface.
      Frame dummyObserver = new Frame();

      //g.setBackground(Color.magenta);
      g.clearRect(0, 0, scaledW, scaledH);
      g.drawImage(bi, 0, 0, scaledW, scaledH, dummyObserver);
      g.dispose();

      return scaledImg;
   } // scaleImage()

//---------------------------------------------------------------------------
   public synchronized static void startClock() {
      time = System.currentTimeMillis();
   }

//---------------------------------------------------------------------------
   public synchronized static long stopClock() {
      return System.currentTimeMillis() - time;
   }

//===========================================================================
   /** 
    * Runnable inner class that loads images
    * 
    * @author Tom Perry
    * @version 
    */
   class ImageLoader extends Observable implements Runnable {
      private String nodeID;
      private Observer requester;

      private ImageLoader(String nodeID, Observer requester) {
	 this.nodeID = nodeID;
	 this.requester = requester;
      } 

      /* 
       * Could return the loaded image in the argument to
       * notifyObservers(), but this would force all receiving
       * classes to be prepared for multithreadedness.
       *
       * Instead, we make the ImageLoader observable and notify
       * requesting classes when the load has finished.  They can
       * then request the loaded image whenever they like
       * (generally "now").
       *
       * We're already getting enough parallelism from loading the
       * images, so a single-threaded getImage() isn't a
       * bottleneck,
       */
      public void run() {
	 addObserver(requester);
	 if (!isLoaded(nodeID)) {
	    startClock();
	    BufferedImage img = fetchImage(nodeID);
	    ImageViewManager.printDebugMessage(
		  "Loaded an image for "+nodeID+" in "+stopClock()+"ms");
	    assert img != null;
	    imageMap.put(nodeID, img);
	    while (!isLoaded(nodeID)) {
	       try { img.wait(10); }
	       catch (InterruptedException ie) {}
	    }
	 } else {
	    ImageViewManager.printDebugMessage("Loading an image for "+nodeID+"... already loaded.");
	 }
	 setChanged();
	 notifyObservers(nodeID);
	 deleteObserver(requester);
      }
   }
}
