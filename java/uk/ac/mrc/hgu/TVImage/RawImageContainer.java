package uk.ac.mrc.hgu.TVImage;

import java.awt.image.*;

/** 
 * Image container subclass that displays raw images.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
public class RawImageContainer extends StandardImageContainer {

//-------------------------------------------------------------------
   public RawImageContainer(String s, BufferedImage i,
	 String c, String emageId, String geneName, int sortIndex,
	 ImageViewManager tvim) {

      super(s,i,c,emageId,geneName,sortIndex,tvim);
   }
//-------------------------------------------------------------------
   public String getHyperlinkTarget() {
      //System.out.println("RawImageContainer: enter getHyperlinkTarget, emageId = "+emageId);
      //String ret = appSettings.getRawImageLocation(this.emageId);
      //System.out.println("RawImageContainer: exit getHyperlinkTarget, location = "+ret);
      return appSettings.getRawImageLocation(this.emageId);
   }
//-------------------------------------------------------------------

}
