//pragma ident "MRC HGU $Id$"
/*!
* \file         WlzObjModel.java
* \author       Nick Burton
* \date         May 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Model class for Wlz objects (MVC)
* \todo         -
* \bug          None known.
*/
package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;

import uk.ac.mrc.hgu.Wlz.*;

//-------------------------------------------------------------
// 'model' class for wlz objects
//-------------------------------------------------------------

/**
 *   <b>Model</b> class for SectionViewer component.
 *   <br>Encapsulates data for 3D Woolz Objects.
 */
public class WlzObjModel implements WlzObjectType {

   /**
    *   File containing a 3D Woolz object.
    */
   File imgFile;

   /**
    *   Stream used to read a 3D Woolz object.
    */
   WlzFileStream in = null;

   /**
    *   3D Woolz object
    */
   WlzObject _obj = null;

   /**
    *   2D section through a 3D Woolz object.
    */
   WlzObject _sectObj = null;
   
   /**
    *   Bounding box of a 3D Woolz object.
    */
   WlzIBox3 _bBox3 = null;


   /**
    *   Flag to control content of debugging output.
    */
   int manyFacts;

//----------------------------------------------
// constructors
//----------------------------------------------
   /**
    *   Constructs a WlzObjModel given a Woolz filename;
    *   @param     imgStr the Woolz object filename.
    */
   public WlzObjModel(String imgStr) {
      manyFacts = 0;
      setWlzObj(imgStr);
   }

   /**
    *   Constructs a WlzObjModel given a Woolz File;
    *   @param     imgFile the Woolz object File.
    */
   public WlzObjModel(File imgFile) {
      manyFacts = 0;
      setWlzObj(imgFile);
   }

//----------------------------------------------
// finalizer
//----------------------------------------------
   /**
    *   Closes WlzFileStream at destruction of the WlzObjModel.
    */
   protected void finalize() {
     try {
	if (in != null) in.close();
     }
     catch (IOException e) {
        System.out.println (e.getMessage());
     }
   }
//----------------------------------------------
   /**
    *   Reads in a 3D Woolz object from the given file.
    *   @param imgFileStr Filename of Woolz object to be read in
    */
   public void setWlzObj(String imgFileStr) {
      try {
	 imgFile = new File(imgFileStr);
	 setWlzObj(imgFile);
      }
      catch (Exception e) {
         System.out.println("setWlzObj");
         System.out.println(e.getMessage());
      }
   }

//----------------------------------------------
   /**
    *   Reads in a 3D Woolz object from the given File.
    *   @param imgFileStr File of Woolz object to be read in
    */
   public void setWlzObj(File imgFile) {

      try {
	 in = new WlzFileInputStream(imgFile.getAbsolutePath());
	 _obj = WlzObject.WlzReadObj(in);
	 _bBox3 = WlzObject.WlzBoundingBox3I(_obj);
      }
      catch (IOException e) {
	System.err.println(e);
      }
      catch (WlzException e) {
	System.err.println("setWlzObj");
	System.err.println(e);
      }
   }

//----------------------------------------------
   /**
    *   Sets _obj to be the given 3D Woolz object.
    *   @param object Woolz object.
    */
   public void setWlzObj(WlzObject object){
     _obj = object;
   }

//----------------------------------------------
   /**
    *   Returns a 2D section through the current grey-level 3D Woolz Object.
    *   @param VS View Structure defining the section.
    *   @return The section.
    */
   public WlzObject makeSection(WlzThreeDViewStruct VS) {

      try {
	 _sectObj = _obj.WlzGetSectionFromObject(
	      _obj,
	      VS,
	      WlzInterpolationType.WLZ_INTERPOLATION_NEAREST);
      }
      catch (WlzException e) {
	System.err.println("makeSection");
	System.err.println(e);
      }

      return _sectObj;

   } // makeSection()

//----------------------------------------------
   /**
    *   Returns a 2D section through the given 3D Woolz Object.
    *   @param VS View Structure defining the section.
    *   @param  obj Woolz object through which to cut a section.
    *   @return The section.
    */
   public WlzObject makeSection(WlzObject obj, WlzThreeDViewStruct VS) {

      WlzObject ret = null;
      try {
	 ret = obj.WlzGetSectionFromObject(
	      obj,
	      VS,
	      WlzInterpolationType.WLZ_INTERPOLATION_NEAREST);

	 if(ret != null) {
	    
	 }
      }
      catch (WlzException e) {
	//System.err.println("makeSection2");
	//System.err.println(e);
      }

      return ret;
   }

//----------------------------------------------
   /**
    *   Returns a thresholded region of a 2D Woolz object.
    *   @param  obj Woolz object to be thresholded.
    *   @param  val threshold (grey) value.
    *   @param  highLow flag which determines if thresholding is to be
    *   above or below <em>val</em>.
    *   @return The result of the threshold operation.
    */
   public WlzObject threshold(WlzObject obj, int val, int highLow) {

      //System.out.println("threshold value = "+val);
      WlzPixelV WPV = new WlzPixelV(val);
      WlzObject ret = null;
      try {
	 ret = obj.WlzThreshold(obj, WPV, highLow);
      }
      catch (WlzException e) {
	System.err.println("threshold");
	System.err.println(e);
      }

      return ret;
   }

//----------------------------------------------
   /**
    *   Returns a smoothed version of a Woolz object.
    *   Used in thresholding operations.
    *   @param  obj Woolz object to be smoothed.
    *   @param wx
    *   @param wy
    *   @param x_deriv
    *   @param y_deriv
    *   @return the smoothed Woolz object.
    */
   public WlzObject smooth(WlzObject obj,
                           double wx,
			   double wy,
			   int x_deriv,
			   int y_deriv) {

      WlzObject ret = null;
      try {
	 ret = obj.WlzGauss2(obj,wx,wy,x_deriv,y_deriv);
      }
      catch (WlzException e) {
	System.err.println("smooth");
	System.err.println(e);
      }

      return ret;
   }

//----------------------------------------------
   /**
    *   Returns the region of a 2D Woolz object which intersects
    *   a given 2D constraint.
    *   @param obj the Woolz object to be constrained.
    *   @param constraint threshold constraint region.
    *   @return the intersecting region.
    */
   public WlzObject constrain(WlzObject obj,
                           WlzObject constraint) {

      if(obj == null) {
	 System.err.println("constrain ... obj == null");
         return null;
      }

      if(constraint == null) {
	 System.err.println("constrain ... constraint == null");
         return obj;
      }

      WlzIBox2 _bBox1 = null; // grey obj
      WlzIBox2 _bBox2 = null; // constraint obj
      try {
	 _bBox1 = WlzObject.WlzBoundingBox2I(obj);
	 _bBox2 = WlzObject.WlzBoundingBox2I(constraint);
      }
      catch(WlzException we) {
	 System.err.print("constrain");
	 System.err.print(we.getMessage());
      }
      /*
      printBoundingBox2(_bBox1);
      printBoundingBox2(_bBox2);
      */

      WlzObject ret = null;
      WlzObject obj1 = null;
      try {
	 obj1 = WlzObject.WlzIntersect2(obj,constraint);
	 /* WLZ_EMPTY_OBJ = 127 */
	 if(obj1 != null) {
	    ret = WlzObject.WlzBndSetGreyValues(obj1, obj);
	 } else {
	    System.err.println("constrain ... obj1 == null");
	 }
      }
      catch (WlzException e) {
	System.err.println("constrain");
	System.err.println(e);
      }

      return ret;
   }

//----------------------------------------------
   /**
    *   Returns the Woolz object from the given array 
    *   that contains the given x,y point.
    *   @param objArray Collection of 2D Woolz objects in a plane.
    *   @param x x coord of test point.
    *   @param y y coord of test point.
    *   @return the containing object.
    */
   private WlzObject contains(WlzObject[] objArray,
				     double x,
		                     double y) {

      WlzObject ret = null;

      int num = objArray.length;
      int isInside = 0;
      int i = -1;

      if(num > 0) {
	 for(i=0; i<num; i++) {
	    try {
	       isInside = WlzObject.WlzInsideDomain(
	                        objArray[i], 0.0, y, x);
	    }
	    catch(WlzException e) {
	       System.out.println("contains "+e.getMessage());
	    }
	    if(isInside != 0) { // assume this means true
	       ret = objArray[i];
	       break;
	    }
	 }
      }

      return ret;

   }

//----------------------------------------------
   /**
    *   takes a 2D Woolz object, decomposes it and returns the region
    *   containing the given point(x,y). 
    *   @param obj the composite Woolz object.
    *   @param x x coord of contained point
    *   @param y y coord of contained point
    *   @return the decomposed region that contains point(x,y).
    */
   public WlzObject select(WlzObject obj,
                           double x,
		           double y) {

      WlzObject ret = null;
      try {
         ret = obj.Wlz2DContains(obj, x, y);
      }
      catch (WlzException e) {
	System.err.println("select");
	 System.err.println(e);
      }

      /*
      WlzObject objArr[][] = new WlzObject[1][];
      int siz[] = new int[1];
      //int maxnum = 1024;
      int maxnum = 2048;
      int ign = 0;
      int connect = 1; // 1=>8-connected, 2=>4-connected

      try {
	 WlzObject.WlzLabel(obj,
			    siz,
			    objArr,
			    maxnum,
			    ign,
			    connect);
      }
      catch (WlzException e2) {
         System.err.println("label");
         System.err.println(e2.getMessage());
      }

      System.err.println("x = "+x+", y = "+y);
      if(objArr[0] == null) {
	System.err.println("objArr[0] == null");
      } else {
	 ret = contains(objArr[0], x, y);
      }
      */

      return ret;
   }

//----------------------------------------------
   /**
    *   Calculates the intersection of a section with
    *   the bounding box specified by the given View Structure.
    *   <em><b>Obsolete</b></em>.
    *   @param viewStr the current View Structure.
    *   @param dstSizeArrayVtxs size of intersection point array.
    *   @param dstArrayVtxs array containing the intersection points.
    *   @return the number of intersection points. Should be 12 unless 
    *   an error has occurred.
    */
   public int Wlz3DViewGetBoundingBoxIntersectionA (
                   WlzThreeDViewStruct viewStr,
                   int [] dstSizeArrayVtxs,
                   WlzDVertex3 [][] dstArrayVtxs) {
      int ret = 0;
      try {
        ret = _obj.Wlz3DViewGetBoundingBoxIntersectionA(
                        viewStr, dstSizeArrayVtxs, dstArrayVtxs);
      } catch (WlzException e) {
	System.err.println("Wlz3DViewGetBoundingBoxIntersectionA");
        System.out.println(e);
      }
      return ret;
   }
//----------------------------------------------
   /**
    *   Returns the current grey-level section.
    *   @return the current grey-level section.
    */
   public WlzObject getSection() {
      return _sectObj;
   }

//----------------------------------------------
   /**
    *   Returns the current grey-level 3D Woolz object.
    *   @return the current grey-level 3D Woolz object.
    */
   public WlzObject getThreeDObj() {
      return _obj;
   }

//----------------------------------------------
   /**
    *   Returns the bounding box for the current grey-level 3D Woolz object.
    *   @return the bounding box for the current grey-level 3D Woolz object.
    */
   public WlzIBox3 getBBox() {
      return _bBox3;
   }

//----------------------------------------------
   /**
    *   Returns the point in 3D Woolz object space
    *   corresponding to a point in the 2D space of the current section.
    *   @param pt the 2D point.
    *   @param VS the current View Structure.
    *   @return the 3D point.
    */
   public double[] get3DPoint(Point pt, WlzThreeDViewStruct VS) {

      double ret[] = new double[3];
      double dist[] = new double[1];

      ret[0] = -1.0;
      ret[1] = -1.0;
      ret[2] = -1.0;

      try {
	 _obj.Wlz3DViewGetDist(VS, dist);
      }
      catch (WlzException e) {
        System.err.println("get3DPoint");
        System.out.println(e.getMessage());
      }
      WlzDVertex3 xy = new WlzDVertex3(pt.getX(), pt.getY(), dist[0]);
      WlzDVertex3 xyz[] = new WlzDVertex3[1];

      try {
	 _obj.Wlz3DSectionTransformInvVtxR(VS, xy, xyz);
      }
      catch (WlzException e) {
	System.err.println("get3DPoint");
         System.out.println(e.getMessage());
      }

      ret[0] = xyz[0].vtX;
      ret[1] = xyz[0].vtY;
      ret[2] = xyz[0].vtZ;
      //...................................
      return ret;
   }
//----------------------------------------------
   /**
    *   Returns the point in the 2D space of the current section
    *   corresponding to a point in the 3D Woolz object space.
    *   @param pt3D the 3D point.
    *   @param VS the current View Structure.
    *   @return the 2D point.
    *   Note that this is a 3 element array, the x & y coordinates 
    *   corresponding to elements 0 & 1 respectively.
    */
   public double[] get2DPoint(double pt3D[], WlzThreeDViewStruct VS) {

      double ret[] = new double[3];
      ret[0] = -1.0;
      ret[1] = -1.0;
      ret[2] = -1.0;
      WlzDVertex3 xyz3D = new WlzDVertex3(pt3D[0],
                                          pt3D[1],
                                          pt3D[2]);
      WlzDVertex3 xyz2D[] = new WlzDVertex3[1];

      try {
         _obj.Wlz3DSectionTransformVtxR(VS, xyz3D, xyz2D);
      }
      catch (WlzException e) {
         System.out.println(e.getMessage());
      }

      ret[0] = xyz2D[0].vtX;
      ret[1] = xyz2D[0].vtY;
      ret[2] = xyz2D[0].vtZ;
      //...................................
      return ret;
   }

//----------------------------------------------
   /**
    *   Returns the max and min x,y,z coordinates of
    *   the bounding box for the current 3D Woolz object.
    *   @param VS the current View Structure.
    *   @return the max & min x,y,z values in a 6 element array.
    *   Elements 0,1,2 correspond to min x,y,z respectively.
    *   Elements 3,4,5 correspond to max x,y,z respectively.
    *   
    */
   public double[] getMaxMin(WlzThreeDViewStruct VS) {

      double ret[] = new double[6];
      double X[] = new double[1];
      double Y[] = new double[1];
      double Z[] = new double[1];
      X[0] = -1.0;
      Y[0] = -1.0;
      Z[0] = -1.0;

      try {
	 _obj.Wlz3DViewGetMaxvals(VS, X, Y, Z);
         ret[0] = X[0];
         ret[1] = Y[0];
         ret[2] = Z[0];
      }
      catch (WlzException e) {
	System.err.println("getMaxMin");
         System.out.println(e.getMessage());
      }

      try {
	 _obj.Wlz3DViewGetMinvals(VS, X, Y, Z);
         ret[3] = X[0];
         ret[4] = Y[0];
         ret[5] = Z[0];
      }
      catch (WlzException e) {
	System.err.println("getMaxMin2");
         System.out.println(e.getMessage());
      }

      //...................................
      return ret;
   }
//----------------------------------------------
   /**
    *   Displays pertinent details
    *   of the current grey-level 3D Woolz object.
    */
   public void printFacts() {

      System.out.println("WLZ FACTS");
      String facts[] = {""} ;
      try {
         WlzObject.WlzObjectFacts(_obj, null, facts, manyFacts);
         System.out.println(facts[0]);
	 System.out.println("*******************************");
	 /*
         WlzObject.WlzObjectFacts(_sectObj, null, facts, manyFacts);
         System.out.println(facts[0]);
	 System.out.println("*******************************");
	 */
      }
      catch (WlzException e) {
	System.err.println("printFacts");
	System.err.println(e);
      }
   } // printFacts()

//----------------------------------------------
   /**
    *   Displays the coordinates of the
    *   bounding box for the current 3D grey-level object.
    */
   public void printBoundingBox() {
      printBoundingBox3(_bBox3);
   }

//----------------------------------------------
   /**
    *   Displays the coordinates of the
    *   bounding box for the given 3D Woolz object.
    *   @param box the bounding box whose coordinates are
    *   to be displayed.
    */
   public void printBoundingBox3(WlzIBox3 box) {
      System.out.println("BOUNDING BOX");
      try {
	 System.out.println("xmin: "+box.xMin+", xmax: "+box.xMax);
	 System.out.println("ymin: "+box.yMin+", ymax: "+box.yMax);
	 System.out.println("zmin: "+box.zMin+", zmax: "+box.zMax);
	 System.out.print("Bounding Box centre, ");
	 System.out.println((box.xMin+(box.xMax-box.xMin)/2.0)+", "+
	                    (box.yMin+(box.yMax-box.yMin)/2.0)+", "+
	                    (box.zMin+(box.zMax-box.zMin)/2.0));
      }
      catch (Exception e) {
	System.err.println(e);
      }
      System.out.println("hhhhhhhhhhhhhhhhhhhhhhhhhhhhhh");
   } // printBoundingBox()

//----------------------------------------------
   /**
    *   Displays the coordinates of the
    *   bounding box for the given 2D Woolz object.
    *   @param box the bounding box whose coordinates are
    *   to be displayed.
    */
   public void printBoundingBox2(WlzIBox2 box) {
      System.out.println("BOUNDING BOX");
      try {
	 System.out.println("xmin: "+box.xMin+", xmax: "+box.xMax);
	 System.out.println("ymin: "+box.yMin+", ymax: "+box.yMax);
	 System.out.print("Bounding Box centre, ");
	 System.out.println((box.xMin+(box.xMax-box.xMin)/2.0)+", "+
	                    (box.yMin+(box.yMax-box.yMin)/2.0));
      }
      catch (Exception e) {
	System.err.println(e);
      }
      System.out.println("hhhhhhhhhhhhhhhhhhhhhhhhhhhhhh");
   } // printBoundingBox()

//-------------------------------------------------------------
// handle all objects that are interested in changes
//-------------------------------------------------------------

   /**
    *   A list of ChangeListeners which are
    *   listening for events fired from the WlzObjModel.
    */
   protected EventListenerList changeListeners =
      new EventListenerList();

   /**
    *   Adds a ChangeListener to the EventListenerList.
    */
   public void addChangeListener(ChangeListener x) {
      changeListeners.add (ChangeListener.class, x);

      // bring it up to date with current state
      x.stateChanged(new ChangeEvent(this));
   }

   /**
    *   Removes a ChangeListener from the EventListenerList.
    */
   public void removeChangeListener(ChangeListener x) {
      changeListeners.remove (ChangeListener.class, x);
   }

   /**   An event that will be fired from WlzObjModel */
   private ChangeEvent ce;
   /**  A local copy of the list of ChangeListeners */
   private Object[] listeners;
   /**  One of the list of Change Listeners */
   private ChangeListener cl;

   /**
    *   Fires a ChangeEvent from the WlzObjModel.
    */
   protected void fireChange() {
      // Create the event:
      ce = new ChangeEvent(this);
      // Get the listener list
      listeners = changeListeners.getListenerList();
      // Process the listeners last to first
      // List is in pairs, Class and instance
      for (int i = listeners.length-2; i >= 0; i -= 2) {
	 if (listeners[i] == ChangeListener.class) {
	    ChangeListener cl = (ChangeListener)listeners[i+1];
	    cl.stateChanged(ce);
	 }
      }
   } // fireChange

//-------------------------------------------------------------
} // class WlzObjModel
