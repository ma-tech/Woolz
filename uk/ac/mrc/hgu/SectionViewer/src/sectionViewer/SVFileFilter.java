package sectionViewer;
import sectionViewer.*;

import java.io.File;
import javax.swing.filechooser.*;

/**
 *   File filter which accepts .jpg or .jpeg files
 */
public class SVFileFilter extends FileFilter {

   /**
    *   Default constructor
    */
   public SVFileFilter() {
   }

   /**
    *   Tests a File type
    *   @param fil the File to be tested
    *   @return true if the file is a .jpg .jpeg or a directory.
    */
   public boolean accept(File fil) {

      boolean ret = false;

      if(fil.isDirectory()) return true;

      String filstr = fil.getName();
      int indx = filstr.lastIndexOf(".");

      if(indx != -1) {
	 String extn = filstr.substring(indx+1);
	 //System.out.println("extn = "+extn);

	 if(extn.equals("jpg")) ret = true;
	 if(extn.equals("jpeg")) ret = true;

      }
      return ret;
   }

   /**
    *   Returns a String description of acceptable File type(s).
    *   @return String description of acceptable File type(s).
    */
   public String getDescription() {
      return "jpg, jpeg";
   }

}



