/* file name  : HeatmapImageFactory.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
 * created    : Mon 03 Sep 2007 09:38:05 BST
 */
package uk.ac.mrc.hgu.TVImage;

import java.awt.image.*;
import java.util.*;

/** 
 * Image factory subclass that produces heatmap images.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
class HeatmapImageFactory extends ImageServer {
   public HeatmapImageFactory(ImageViewManager tvim, XmlAppSettings xml) {
      super(tvim, xml);
   }

   public String getImagePath(String nodeID) {
      return xml.getHeatmapThumbnailLocation(nodeID);
   }

   public HeatmapImageContainer getMedia(String id) {
      if (tvim.isDebugMode()) {
	 startClock();
      }

      int sortIndex = tvim.getInorderSearchIndex(id);
      BufferedImage img = getImage(id);
      //System.out.println("HeatmapImageFactory: getMedia id = "+id);
      //System.out.println("buffered image size = "+img.getWidth()+", x "+img.getHeight());

      if (tvim.isDebugMode()) {
	 ImageViewManager.printDebugMessage(
	       "Got heatmap image for "+id+" in "+stopClock()+"ms");
      }

      return new HeatmapImageContainer(id, img, tvim.getCorrelation(id),
	    sortIndex, xml.getHeatmapVideoLocation(id));
   }
}
