/**
 * 
 */

package uk.ac.mrc.hgu.TVImage;

import javax.xml.parsers.*;

import org.xml.sax.*;
import org.xml.sax.helpers.*;
import org.w3c.dom.*;

import java.io.*;
import java.util.*;

/**
 * @author paul smith
 *
 */
public class XmlAppSettings {
	private HashMap<String, String> appSettings;
	private boolean parseSettingsForThisApp;
	
	
	/**
	 * The constructor will take the filename and parse it using
	 * the DOM framework.  Until the supplied appName is found as
	 * an attribute of an ELEMENT named App nothing will be stored
	 * All document elements under the App tag will be stored in a
	 * HashMap (tagname as the key, TEXT node under it as value)
	 * Once moved off the relevant App tag don't store any further
	 * values
	 * @param filename
	 * @param appName
	 */
	public XmlAppSettings(String filename, String appName) throws Exception {

	   //System.out.println("XmlAppSettings: enter constructor filename = "+filename+", appName = "+appName);

	   //initialise storage of settings
	   appSettings = new HashMap<String, String>(); 

	   parseSettingsForThisApp = false;

	   //Step 1: create a DocumentBuilderFactory
	   DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();


	   //Step 2: Create a DocumentBuilder
	   DocumentBuilder db = null;
	   try
	   {
	      db = dbf.newDocumentBuilder();
	   }
	   catch (ParserConfigurationException pce)
	   {
	      System.err.println(pce);
	      throw pce;
	   }

	   //Step 3: parse the input file
	   Document doc = null;
	   if (filename == null)
	      throw new IllegalArgumentException();

	   try
	   {
	      if (filename.startsWith("http"))
		 doc = db.parse(filename);
	      else
		 doc = db.parse(new File(filename));
	   }
	   catch (SAXException se)
	   {
	      System.err.println("Could not parse XML file: "+se.getMessage());
	      throw se;
	   }

	   //Step 4: echo the document
	   if (doc != null)
	      echo(doc, appName);

	   //printOut();
	   //System.out.println("XmlAppSettings: exit constructor");
	}
	
	/**
	 * this will walk the tree recursively.
	 * Once it finds the specified application name it sets a boolean to denote that it has found the desired app
	 * Then it stores the other elements and their children using storeElementAndValueInHashMap
	 * Once it finds a different application it stops storing the elements
	 * @param n node under examination
	 * @param appName name of the applcation for which we are looking for values
	 */
	private void echo(Node n, String appName)
	{
		if (n instanceof Element) {
			if (n.getNodeName().equals("Apps")) {
				//do nothing
			}
			else if (n.getNodeName().equals("App")) {
				//want to grab the app name
				//if it matches then want to grab the rest of the settings
				parseSettingsForThisApp = true;
				//otherwise ignore it
			}
			else if (parseSettingsForThisApp) {
				//want to store the element and its value in the map
				for (Node child= n.getFirstChild(); child != null; child = child.getNextSibling()) {
					storeElementAndValueInHashMap(n.getNodeName(), child);
				}
			}
		}
		for (Node child = n.getFirstChild(); child != null; child = child.getNextSibling()) {
			echo(child, appName);
		}
	}
	
	/**
	 * If childNode is a Text element then store it as a value in the hash table using name of parent node as index 
	 * @param elementName name of parent node
	 * @param childNode
	 */
	private void storeElementAndValueInHashMap(String elementName, Node childNode)
	{
		if (childNode instanceof Text) {
			appSettings.put(elementName, childNode.getNodeValue());
		}
	}
	
	/**
	 * for debugging
	 * will walk the hash map and print out all the values
	 *
	 */
	public void printOut()
	{
		Iterator it = appSettings.entrySet().iterator();
		while (it.hasNext()) {
	        Map.Entry pairs = (Map.Entry)it.next();
	        System.out.println(pairs.getKey() + " = " + pairs.getValue());
	    }
	}
	
	
	/**
	 * 
	 * @return the location for the atr file
	 */
	public String getATRFileLocation()
	{
		return appSettings.get("ATRFileLocation");
	}
	
	/**
	 * 
	 * @return the location for the cdt file
	 */
	public String getCDTFileLocation()
	{
		return appSettings.get("CDTFileLocation");
	}
	
	/**
	 * 
	 * @return the location for the gtr file
	 */
	public String getGTRFileLocation()
	{
		return appSettings.get("GTRFileLocation");
	}
	
	/**
	 * 
	 * @return the location for the jtv file
	 */
	public String getJTVFileLocation()
	{
		return appSettings.get("JTVFileLocation");
	}
	
	/**
	 * 
	 * @return the placeholder where the unique id for the images is to be substituted into the file locations
	 */
	public String getPlaceholderForSpecificImage()
	{
		return appSettings.get("PlaceholderForSpecificImage");
	}
	
	/**
	 * 
	 * @param specificId emage id for this gene
	 * @return url for the emage page for the supplied emage id
	 */
	public String getEmagePageLocation(String specificId)
	{
		return substituteValueForPlaceholder("EmagePageURL", specificId);
	}
	
	/**
	 * 
	 * @param geneName name of the gene
	 * @return url for the search for the supplied gene name
	 */
	public String getGeneNameQueryString(String geneName)
	{
		return substituteValueForPlaceholder("GeneNameQueryURL", geneName);
	}
	
	/**
	 * 
	 * @param specificImageId unique identifier for the specified image
	 * @return location for raw thumbnail for specified id
	 */
	public String getRawThumbnailLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("RawThumbnailURL", specificImageId);
	}
	
	/**
	 * 
	 * @param specificImageId specificImageId unique identifier for the specified image
	 * @return location for raw image for specified id
	 */
	public String getRawImageLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("RawImageURL", specificImageId);
	}
	
	/**
	 * 
	 * @param specificImageId ImageId unique identifier for the specified image
	 * @return location for annotated thumbnail for specified id
	 */
	public String getAnnotatedThumbnailLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("AnnotatedThumbnailURL", specificImageId);
	}
	
	/**
	 * 
	 * @param specificImageId unique identifier for the specified image
	 * @return location for annotated image for specified id
	 */
	public String getAnnotatedImageLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("AnnotatedImageURL", specificImageId);
	}
	
	/**
	 * 
	 * @param specificImageId unique identifier for the specified image
	 * @return location for heatmap thumbnail for specified id
	 */
	public String getHeatmapThumbnailLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("HeatmapThumbnailURL", specificImageId);
	}
	
	/**
	 * 
	 * @param specificImageId unique identifier for the specified image
	 * @return location for heatmap image for specified id
	 */
	public String getHeatmapImageLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("HeatmapImageURL", specificImageId);
	}
		
	/**
	 * 
	 * @param specificImageId unique identifier for the specified image
	 * @return thumbnail location for heatmap video for specified id
	 */
	public String getHeatmapVideoThumbnailLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("HeatmapVideoThumbnailURL", specificImageId);
	}

	/**
	 * 
	 * @param specificImageId unique identifier for the specified image
	 * @return Heatmap video location for specified id
	 */
	public String getHeatmapVideoLocation(String specificImageId)
	{
		return substituteValueForPlaceholder("HeatmapVideoURL", specificImageId);
	}	
	
	/**
	 * 
	 * @param keyValue key to string value we want to update
	 * @param replacementValue value placeholder is to be replaced by
	 * @return value stored in hashmap under specified key but with the placeholder replaced with the supplied replacementValue
	 */
	private String substituteValueForPlaceholder(String keyValue, String replacementValue) {
	        //System.out.println("XmlAppSettings: enter substituteValueForPlaceholder");
		String returnValue;
		//System.out.println("Key is " + keyValue);
		returnValue = appSettings.get(keyValue);
		//System.out.println("Value was " + returnValue + " and replacementValue " + replacementValue);
		if (returnValue != null)
			returnValue = returnValue.replace(getPlaceholderForSpecificImage(), replacementValue);
		//System.out.println("return value is " + returnValue);
	        //System.out.println("XmlAppSettings: exit substituteValueForPlaceholder");
		return returnValue;
	}
}
