package sectionViewer;
import sectionViewer.*;

import java.io.File;
import javax.swing.filechooser.*;

public class SVFileFilter extends FileFilter {


   public SVFileFilter() {
   }

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

   public String getDescription() {

      return "jpg, jpeg";
   }

}



