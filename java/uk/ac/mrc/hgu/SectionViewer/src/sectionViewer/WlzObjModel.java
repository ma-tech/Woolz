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
public class WlzObjModel implements WlzObjectType {

   File imgFile;

   WlzFileStream in = null;
   WlzObject _obj = null;
   WlzObject _sectObj = null;
   WlzIBox3 _bBox3 = null;

   int manyFacts;

//----------------------------------------------
// constructors
//----------------------------------------------
   public WlzObjModel(String imgStr) {
      manyFacts = 0;
      setWlzObj(imgStr);
   }

   public WlzObjModel(File imgFile) {
      manyFacts = 0;
      setWlzObj(imgFile);
   }

//----------------------------------------------
// finalizer
//----------------------------------------------
   protected void finalize() {
     try {
	if (in != null) in.close();
     }
     catch (IOException e) {
        System.out.println (e.getMessage());
     }
   }
//----------------------------------------------
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

   public void setWlzObj(WlzObject object){
     this._obj = object;
   }
//----------------------------------------------
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
   public WlzObject getSection() {
      return _sectObj;
   }

//----------------------------------------------
   public WlzObject getThreeDObj() {
      return _obj;
   }

//----------------------------------------------
   public WlzIBox3 getBBox() {
      return _bBox3;
   }

//----------------------------------------------
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
   public void printBoundingBox() {
      printBoundingBox3(_bBox3);
   }

//----------------------------------------------
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

   // keep track of all the listeners to this model
   protected EventListenerList changeListeners =
      new EventListenerList();

   // add a listener to the register
   public void addChangeListener(ChangeListener x) {
      changeListeners.add (ChangeListener.class, x);

      // bring it up to date with current state
      x.stateChanged(new ChangeEvent(this));
   }

   // remove a listener from the register
   public void removeChangeListener(ChangeListener x) {
      changeListeners.remove (ChangeListener.class, x);
   }

   private ChangeEvent ce;
   private Object[] listeners;
   private ChangeListener cl;

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
} // class ViewStruct
