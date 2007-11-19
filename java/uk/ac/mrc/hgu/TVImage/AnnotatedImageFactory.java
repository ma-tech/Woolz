/* file name  : AnnotatedImageFactory.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
 * created    : Mon 03 Sep 2007 09:12:24 BST
 */
package uk.ac.mrc.hgu.TVImage;

import java.awt.image.*;
import java.util.*;

/** 
 * Image factory subclass that produces annotated images.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
class AnnotatedImageFactory extends ImageServer {

   public AnnotatedImageFactory(ImageViewManager tvim, XmlAppSettings xml) {
      super(tvim, xml);
   }

   protected String getImagePath(String nodeID) {
      String emageId = tvim.getEMAGEID(nodeID);
      //strip off the substring EMAGE: that starts the EmageID
      if (emageId.startsWith("EMAGE:") && emageId.length() > 6) {
	 emageId = emageId.substring(6);
      }
      return xml.getAnnotatedThumbnailLocation(emageId);
   }

   public AnnotatedImageContainer getMedia(String id) {
      //System.out.println("AnnotatedImageFactory: enter getMedia "+id);
      ImageViewManager.printDebugMessage("Getting annotated image for "+id);

      String emageId = tvim.getEMAGEID(id);
      String geneName = tvim.getGeneName(id);
      BufferedImage img = getImage(id);

      //System.out.println("AnnotatedImageFactory: exit getMedia "+id);
      return new AnnotatedImageContainer(id, img,tvim.getCorrelation(id),
	    emageId, geneName, tvim.getGeneIndex(id), tvim);
   }
}
