package sectionViewer;
import sectionViewer.*;

import java.net.*;
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.help.*;

/*
   based on comp.lang.java.programmer newsgroup reply from
   Richard Trahan (rtrahan@monmouth.com) 2001-05-01
*/
public class HelpUtils {

    /** Returns the full path name of the directory from which the
     * class of the argument object was loaded. Makes use of the fact
     * that the class loader stores that directory into the reflection
     * data of the class. Use this when you want to open files which
     * cannot have package names, like .bmp, .gif, .txt, etc. Be sure
     * this class is in the same directory as those files.
     * @return The full path, like d:\com\rockyboat\xutils. Note
     * that the correct platform-dependent file separator will be used
     * (\ on Win32, / on unix, etc.), and that the path does not end 
     * with the separator character.
     * @param obj The calling class whose full path origin is to be
     * determined.
     */
    public static String getLoadPath(Object obj) {
        Class c = obj.getClass(); // gives us access to the class loader's data
        String s = c.getName(); // full package name
        s = s.substring(s.lastIndexOf('.')+1); // works even when lastIndexOf returns -1
            // s = name of the calling class, without package (like "myclass")
            // (we had to remove the package path because getResource puts it back)
        URL u = c.getResource(s+".class"); // like myclass.class
            // getResource is guaranteed to access the directory from which
            // the calling class was loaded
        s = u.getPath();
            // always uses "/", which is a url separator, unrelated to the platform
        int iLast = s.lastIndexOf('/');
        s = s.substring(0,iLast);
            // i.e., the full path to the directory where the calling class
            // was loaded
	/*
	System.out.println("getLoadPath returning "+s);
	if package is in a jar file, returns something like
	file:/export/data0/nickb/MouseAtlas/MyBuild/test/my.jar!/MyPackage
	otherwise
	/export/data0/nickb/MouseAtlas/MyBuild/test/MyPackage
	*/
        return s;
    }

    /*
     Assumes that the help system is in the package directory 
     of the class that is calling it
     ie. calling class is in  package MyPackage
         help files are in directory MyPackage/helpDir
    */
    public static URL getHelpSetURL(Object obj,
				    String helpsetPath)
		            throws MalformedURLException {
        URL hsURL = null;
        String path = getLoadPath(obj);
        if(path.indexOf(".jar!") != -1) { // we're running from a jar
            hsURL = new URL("jar:"+path+"/"+helpsetPath);
        }
        else { 
            hsURL = obj.getClass().getResource(helpsetPath);
        }
	/*
	System.out.println("getHelpSetURL returning "+hsURL);
	if package is in a jar file, returns something like
	jar:file:/export/data0/nickb/MouseAtlas/MyBuild/test/my.jar!/MyPackage/help/ViewerHelp.hs
	otherwise returns something like
	file:/export/data0/nickb/MouseAtlas/MyBuild/test/MyPackage/help/ViewerHelp.hs
	*/
        return hsURL;
    }
}
