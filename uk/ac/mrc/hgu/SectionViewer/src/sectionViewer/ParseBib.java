package sectionViewer;

import java.io.*;
import java.util.*;

/**
 *   Reads View Structure information from a previously saved .bib file.
 */
public class ParseBib {


   //------------------------------------------------
   SectionViewer _SV = null;
   BufferedReader reader = null;
   String line = "";
   String params = "";
   boolean _OK;
   float _fval;
   float _fval1;
   float _fval2;
   float _fval3;
   //------------------------------------------------

   /**
    *   Constructs a ParseBib object for the given SectionViewer.
    *   @param SV the SectionViewer that will use the information
    *   read in from the .bib file.
    */
   protected ParseBib(SectionViewer SV) {
      _OK = false;
      _fval = 0.0f;
      _fval1 = 0.0f;
      _fval2 = 0.0f;
      _fval3 = 0.0f;
      _SV = SV;
   }

   //------------------------------------------------
   /**
    *   Wrapper for bibFileOK(File filename).
    *   @param filename the filename to read.
    */
   protected boolean bibFileOK(String filename) throws Exception {
      return bibFileOK(new File(filename));
   }

   //------------------------------------------------
   /**
    *   Checks the bib file.
    *   @return true if there is an @Wlz3DSectionViewParams section
    */
   protected boolean bibFileOK(File bibFile) throws Exception {

      if(reader == null) {
	 reader = new BufferedReader(new FileReader(bibFile));
	 try {

	    while (reader.ready()) {
	       line = reader.readLine().trim();

	       if (line.startsWith("@Wlz3DSectionViewParams")){
		  _OK = true;
		  break;
	       }
	    }
	 } // try
	 catch (Exception exp) {
	    exp.printStackTrace();
	 }
      }

      return _OK;

   } // bibFileOK()

   //------------------------------------------------
   /**
    *   Wrapper for doParse(File filename).
    *   @param filename the filename to read.
    */
   protected void doParse(String filename) throws Exception {
      doParse(new File(filename));
   }

   //------------------------------------------------
   /**
    *   Parses the bib file.
    *   @param filename the File to read.
    */
   protected void doParse(File bibFile) throws Exception {

      if(!_OK) {
         if(!bibFileOK(bibFile)) return;
      }

      try {
         if(reader == null) {
	    reader = new BufferedReader(new FileReader(bibFile));
	 }

	 StringTokenizer tokeniser = null;

	 while (reader.ready()) {
	    line = reader.readLine().trim();

	    if (line.startsWith("FixedPoint")){
	       params = line.substring(line.indexOf("{") + 1, line.indexOf("}"));
	       tokeniser = new StringTokenizer(params.trim(), ",");
               _fval1 = Float.parseFloat(tokeniser.nextToken());
               _fval2 = Float.parseFloat(tokeniser.nextToken());
               _fval3 = Float.parseFloat(tokeniser.nextToken());
	       _SV.getViewStructModel().setFixedPoint(_fval1,_fval2,_fval3);
	    } else if (line.startsWith("Distance")){
	       params = line.substring(line.indexOf("{") + 1, line.indexOf("}"));
               _fval = Float.parseFloat(params);
               _SV.getDistSetter().setValue(_fval);
	    }else if (line.startsWith("Pitch")){
	       params = line.substring(line.indexOf("{") + 1, line.indexOf("}"));
               _fval = Float.parseFloat(params);
               _SV.getPitchSetter().setValue(_fval);
	    }else if (line.startsWith("Yaw")){
	       params = line.substring(line.indexOf("{") + 1, line.indexOf("}"));
               _fval = Float.parseFloat(params);
               _SV.getYawSetter().setValue(_fval);
	    }else if (line.startsWith("Roll")){
	       params = line.substring(line.indexOf("{") + 1, line.indexOf("}"));
               _fval = Float.parseFloat(params);
               _SV.getRollSetter().setValue(_fval);
	    }else if (line.startsWith("Scale")){
	       params = line.substring(line.indexOf("{") + 1, line.indexOf("}"));
               _fval = Float.parseFloat(params);
               _SV.getZoomSetter().setValue((int)(100.0f*_fval));
	    }else if (line.startsWith("ViewMode")){
	       params = line.substring(line.indexOf("{") + 1, line.indexOf("}"));
               _SV.setViewMode(params);
	    }
	 }
      } // try
      catch (Exception exp) {
	 exp.printStackTrace();
      }
   } // doParse()

} // class ParseBib
