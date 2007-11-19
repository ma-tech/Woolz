/* file name  : AnnotatedImageContainer.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
 * created    : Mon 03 Sep 2007 09:08:06 BST
 * copyright  : MRC HGU
 */
package uk.ac.mrc.hgu.TVImage;

import java.awt.image.*;

/** 
 * Image container subclass that displays annotated images.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 */
public class AnnotatedImageContainer extends StandardImageContainer {

   //-------------------------------------------------------------------
   public AnnotatedImageContainer(String s, BufferedImage i,
	 String c, String emageId, String geneName, int sortIndex,
	 ImageViewManager tvim) {

      super(s,i,c,emageId,geneName,sortIndex,tvim);
   }
   //-------------------------------------------------------------------
   public String getHyperlinkTarget() {
      return appSettings.getAnnotatedImageLocation(this.emageId);
   }
   //-------------------------------------------------------------------
}
