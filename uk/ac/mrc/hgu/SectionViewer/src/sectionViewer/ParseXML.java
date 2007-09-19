package sectionViewer;

// JAXP packages
import java.io.*;
import javax.xml.parsers.*;

import org.xml.sax.*;
import org.xml.sax.helpers.*;


/**
 *   Reads View Structure information from a previously saved .xml file.
 */
public class ParseXML extends DefaultHandler {

    boolean OK;
    String _sval;
    float _fval;
    float _fval1;
    float _fval2;
    float _fval3;
    //ViewStructModel _VSM = null;
    SectionViewer _SV = null;

    //public ParseXML(ViewStructModel VSM) {
    /**
     *   Constructs a ParseXML object for the given SectionViewer.
     *   @param SV the SectionViewer that will use the information
     *   read in from the .xml file.
     */
    public ParseXML(SectionViewer SV) {
       OK = false;
       _sval = "";
       _fval = 0.0f;
       _fval1 = 0.0f;
       _fval2 = 0.0f;
       _fval3 = 0.0f;
       _SV = SV;
    }

    /**
     *   Parser calls this once at the beginning of a document.
     */
    public void startDocument() throws SAXException {
        //System.out.println("starting document");
    }

    /**
     *   Parser calls this for each element in a document.
     */
    public void startElement(String namespaceURI, String localName,
                             String qName, Attributes atts)
        throws SAXException
    {
       //System.out.println("starting element "+localName);
       if(localName.equals("view") == true) {
          OK = true;
          //_SV.getDistSetter().setValue(100.0f); // works !!
       }
    }

    /**
     *   Parser calls this for each element in a document.
     */
    public void characters(char[] buf, int start, int len)
        throws SAXException
    {
       _sval = new String(buf,start,len);
    }

    /**
     *   Parser calls this for each element in a document.
     */
    public void endElement(String namespaceURI, String localName,
                             String qName)
        throws SAXException
    {
       if(OK == true) {
	  if(localName.equals("dist") == true) {
	     _fval = convert(_sval);
	     _SV.getDistSetter().setValue(_fval);
	  } else if(localName.equals("theta") == true) {
	     _fval = convert(_sval);
	     _SV.getYawSetter().setValue(_fval);
	  } else if(localName.equals("phi") == true) {
	     _fval = convert(_sval);
	     _SV.getPitchSetter().setValue(_fval);
	  } else if(localName.equals("zeta") == true) {
	     _fval = convert(_sval);
	     _SV.getRollSetter().setValue(_fval);
	  } else if(localName.equals("mode") == true) {
	     _SV.setViewMode(_sval);
	  } else if(localName.equals("x") == true) {
	     _fval1 = convert(_sval);
	  } else if(localName.equals("y") == true) {
	     _fval2 = convert(_sval);
	  } else if(localName.equals("z") == true) {
	     _fval3 = convert(_sval);
	  } else if(localName.equals("fixed") == true) {
	     _SV.getViewStructModel().setFixedPoint(_fval1,_fval2,_fval3);
	  } else if(localName.equals("viewStruct") == true) {
	     //
	  } else if(localName.equals("zoom") == true) {
	     _SV.getZoomSetter().setValue((Integer.valueOf(_sval)).intValue());
	  } else if(localName.equals("view") == true) {
	     OK = false;
	  } else {
	     System.out.println("unrecognised element in xml file");
	     System.out.println(localName);
	  }
       }
    }

    /**
     *   Parser calls this once after parsing a document.
     */
    public void endDocument() throws SAXException {
        //System.out.println("ending document");
    }

    /**
     *   Converts a numeric string to a float.
     *   @param str the numeric string to convert.
     *   @return a float.
     */
    protected float convert(String str) {
       float ret;
       ret = new Float(str).floatValue();
       return ret;
    }

    /**
     * Convert from a filename to a file URL.
     */
    protected static String convertToFileURL(String filename) {
        // On JDK 1.2 and later, simplify this to:
	String path = "";
	try {
	   path = new File(filename).toURL().toString();
	 }
	 catch(java.net.MalformedURLException e) {
           System.out.println(e.getMessage());
	 }
        return path;
    }

    /**
     * Convert from a filename to a file URL.
     */
    protected static String convertToFileURL(File file) {
	String path = "";
	try {
	   path = file.toURL().toString();
	 }
	 catch(java.net.MalformedURLException e) {
           System.out.println(e.getMessage());
	 }
        return path;
    }

    /**
     *   Wrapper for doParse(File filename).
     */
    public void doParse(String filename) throws Exception {
       doParse(new File(filename));
    }

    /**
     *   Parses the xml file.
     *   @param filename the File to parse.
     */
    public void doParse(File filename) throws Exception {

       //_SV.getDistSetter().setValue(100.0f); // works !!

        // Create a JAXP SAXParserFactory and configure it
        SAXParserFactory spf = SAXParserFactory.newInstance();

        // Set namespaceAware to true to get a parser that corresponds to
        // the default SAX2 namespace feature setting.  This is necessary
        // because the default value from JAXP 1.0 was defined to be false.
        spf.setNamespaceAware(true);

        // Validation part 1: set whether validation is on
        spf.setValidating(false);

        // Create a JAXP SAXParser
        SAXParser saxParser = spf.newSAXParser();

        // Get the encapsulated SAX XMLReader
        XMLReader xmlReader = saxParser.getXMLReader();

        // Set the ContentHandler of the XMLReader
        //xmlReader.setContentHandler(new ParseXML());
        xmlReader.setContentHandler(this);

        // Set an ErrorHandler before parsing
        xmlReader.setErrorHandler(new MyErrorHandler(System.err));

        // Tell the XMLReader to parse the XML document
        xmlReader.parse(convertToFileURL(filename));
    }

    /**
     *   Error handler to report errors and warnings
     */
    private static class MyErrorHandler implements ErrorHandler {
        /** Error handler output goes here */
        private PrintStream out;

        MyErrorHandler(PrintStream out) {
            this.out = out;
        }

        /**
         * Returns a string describing parse exception details
         */
        private String getParseExceptionInfo(SAXParseException spe) {
            String systemId = spe.getSystemId();
            if (systemId == null) {
                systemId = "null";
            }
            String info = "URI=" + systemId +
                " Line=" + spe.getLineNumber() +
                ": " + spe.getMessage();
            return info;
        }

        // The following methods are standard SAX ErrorHandler methods.
        // See SAX documentation for more info.

        public void warning(SAXParseException spe) throws SAXException {
            out.println("Warning: " + getParseExceptionInfo(spe));
        }

        public void error(SAXParseException spe) throws SAXException {
            String message = "Error: " + getParseExceptionInfo(spe);
            throw new SAXException(message);
        }

        public void fatalError(SAXParseException spe) throws SAXException {
            String message = "Fatal Error: " + getParseExceptionInfo(spe);
            throw new SAXException(message);
        }
    }
}



