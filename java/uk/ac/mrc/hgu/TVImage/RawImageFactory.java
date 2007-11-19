package uk.ac.mrc.hgu.TVImage;

import java.awt.image.*;
import java.util.*;

/** 
 * Image factory subclass that produces raw images.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
class RawImageFactory extends ImageServer {

   public RawImageFactory(ImageViewManager tvim, XmlAppSettings xml) {
      super(tvim, xml);
   }

   protected String getImagePath(String nodeID) {
      String emageId = tvim.getEMAGEID(nodeID);
      //strip off the substring EMAGE: that starts the EmageID
      if (emageId.startsWith("EMAGE:") && emageId.length() > 6) {
	 emageId = emageId.substring(6);
      }
      return xml.getRawThumbnailLocation(emageId);
   }

   public RawImageContainer getMedia(String id) {
      //System.out.println("RawImageFactory: enter getMedia "+id);
      ImageViewManager.printDebugMessage("Getting raw image for "+id);

      String emageId = tvim.getEMAGEID(id);
      String geneName = tvim.getGeneName(id);
      BufferedImage img = getImage(id);

      //System.out.println("RawImageFactory: exit getMedia "+id);
      return new RawImageContainer(id, img,tvim.getCorrelation(id),
	    emageId, geneName, tvim.getGeneIndex(id), tvim);
   }
}
