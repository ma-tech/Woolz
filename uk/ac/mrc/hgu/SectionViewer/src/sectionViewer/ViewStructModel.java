//pragma ident "MRC HGU $Id$"
/*!
* \file         ViewStructModel.java
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
// 'model' class for view structure onto wlz objects
//-------------------------------------------------------------

/**
 *   <b>Model</b>  class for SectionViewer.
 *   <br>Encapsulates the View Structure which represents arbitrary
 *   sections through a 3D Woolz object.
 */
public class ViewStructModel implements WlzObjectType {

   /**   Model for the underlying 3D Woolz object */
   WlzObjModel _objModel = null;

   /**   The underlying 3D Woolz object */
   WlzObject _obj = null;

   /**   View Structure for the underlying 3D Woolz object */
   WlzThreeDViewStruct _VS = null;

   /**   The initial fixed point.*/
   double _fpInitial[] = null;
   /**   Normalised vector used in <em>fixed line</em> calculations..*/
   double _norm2[] = null;
   /**   Normalised vector used in <em>fixed line</em> calculations..*/
   double _norm3[] = null;

   /**   Allows pointer values to be retrieved from native C functions. */
   int iValArray[];

   /**   Allows pointer values to be retrieved from native C functions. */
   double dValArray[];
   /**   Allows pointer values to be retrieved from native C functions. */
   double dValArray1[];
   /**   Allows pointer values to be retrieved from native C functions. */
   double dValArray2[];
   /**   Allows pointer values to be retrieved from native C functions. */
   double dValArray3[];

   /**   Toggles display of debugging messages. */
   private final boolean _debug = false;
   /**   True if the 2nd fixed point that defines the fixed line has been defined. */
   private boolean _axisDefined = false;

//----------------------------------------------
// constructors
//----------------------------------------------
   /**
    *   Creates a ViewStructModel for the given WlzObjModel.
    *   @param objModel the given WlzObjModel.
    */
   public ViewStructModel(WlzObjModel objModel) {

      if(_debug) System.out.println("enter ViewStructModel()");
      if(objModel == null) {
         System.out.println("objModel == null");
	 return;
      }
      _objModel = objModel;
      _obj = _objModel.getThreeDObj();
      try {
	 _VS = _obj.WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT);
      }
      catch (WlzException e) {
         System.out.println("ViewStructModel");
         System.out.println(e.getMessage());
      }

      iValArray = new int[1];
      iValArray[0] = 0;
      dValArray = new double[1];
      dValArray1 = new double[1];
      dValArray2 = new double[1];
      dValArray3 = new double[1];
      dValArray[0] = 0.0;
      dValArray1[0] = 0.0;
      dValArray2[0] = 0.0;
      dValArray3[0] = 0.0;

      setInitialViewStruct();
      if(_debug) System.out.println("exit ViewStructModel()");
   }

//----------------------------------------------
//----------------------------------------------
   /**
    *   Initialises the View Structure associated with this ViewStructModel.
    */
   public void setInitialViewStruct() {

      if(_debug) System.out.println("enter setInitialViewStruct()");
      WlzIBox3 bBox3 = _objModel.getBBox();

      _fpInitial = new double[3];

      _fpInitial[0] = bBox3.xMin+(bBox3.xMax-bBox3.xMin)/2.0;
      _fpInitial[1] = bBox3.yMin+(bBox3.yMax-bBox3.yMin)/2.0;
      _fpInitial[2] = bBox3.zMin+(bBox3.zMax-bBox3.zMin)/2.0;

      try {
         _obj.Wlz3DViewSetViewMode(_VS, WlzThreeDViewMode.WLZ_UP_IS_UP_MODE);
         _obj.Wlz3DViewSetUp(_VS, 0.0, 0.0, -1.0);
         _obj.Wlz3DViewSetTheta(_VS, 0.0);
         _obj.Wlz3DViewSetPhi(_VS, 0.0);
         _obj.Wlz3DViewSetDist(_VS,0.0);
         _obj.Wlz3DViewSetFixed(_VS,_fpInitial[0],_fpInitial[1],_fpInitial[2]);
         _obj.Wlz3DViewSetFixed2(_VS,_fpInitial[0],_fpInitial[1],_fpInitial[2]);
         _obj.WlzInit3DViewStruct(_VS,_obj);
         /* adapters not in place yet so don't do fireChange() */
      }
      catch (WlzException e2) {
        System.err.println("setInitialViewStruct");
        System.err.println(e2);
      }
      if(_debug) System.out.println("exit setInitialViewStruct()");

      _norm2 = new double[3];
      _norm3 = new double[3];
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>ViewMode</em> parameter.
    *   @param modeStr a description of the ViewMode.
    *   This must be one of
    *   <ul><li>absolute</li><li>uo is up</li><li>fixed line</li></ul>
    */
   public void setViewMode(String modeStr) {
     int mode;

     if ((modeStr.equals("absolute")) ||
	 (modeStr.equals("Absolute"))) {
       mode = WlzThreeDViewMode.WLZ_ZETA_MODE;
     } else if((modeStr.equals("up is up")) ||
	 (modeStr.equals("Up is up"))) {
       mode = WlzThreeDViewMode.WLZ_UP_IS_UP_MODE;
     } else if((modeStr.equals("fixed line")) ||
	 (modeStr.equals("Fixed line"))) {
       mode = WlzThreeDViewMode.WLZ_FIXED_LINE_MODE;
     } else {
       mode = WlzThreeDViewMode.WLZ_UP_IS_UP_MODE; // default
     }

     try {
       _obj.Wlz3DViewSetViewMode(_VS, mode);
       _obj.WlzInit3DViewStruct(_VS,_obj);
     } catch (WlzException e) {
       System.err.println("setViewMode");
       System.err.println(e);
     }
     _objModel.makeSection(_VS);
     fireChange();
   }

   /**   Conversion from degrees to radians. */
   private double toRads = Math.PI/180.0;
//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Theta</em> parameter.
    *   This corresponds to Yaw or rotation around the Z axis.
    *   @param degs the rotation angle in degrees.
    */
   public void setThetaDeg(double degs) {
      setTheta(degs * toRads);
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Phi</em> parameter.
    *   This corresponds to Pitch or rotation around the X axis.
    *   @param degs the rotation angle in degrees.
    */
   public void setPhiDeg(double degs) {
      setPhi(degs * toRads);
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Zeta</em> parameter.
    *   This corresponds to Roll or rotation around the Y axis.
    *   @param degs the rotation angle in degrees.
    */
   public void setZetaDeg(double degs) {
      setZeta(degs * toRads);
   }
//----------------------------------------------
   /**
    *   Sets the ViewStructModel's <em>Psi</em> parameter.
    *   Note that this is <b>not</b> a View Structure parameter,
    *   it is only used by the ViewStructModel to calculate
    *   the change in Theta Phi and Zeta required to produce
    *   a rotation of Psi around the fixed line..
    *   @param degs the rotation angle in degrees.
    */
   public void setPsiDeg(double degs) {
      setPsi(degs * toRads);
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Theta</em> parameter.
    *   This corresponds to Yaw or rotation around the Z axis.
    *   @param rads the rotation angle in radians.
    */
   public void setTheta(double rads) {
     try {
       _obj.Wlz3DViewSetTheta(_VS, rads);
       _obj.WlzInit3DViewStruct(_VS,_obj);
     }
     catch (WlzException e) {
       System.err.println("setTheta");
       System.err.println(e);
     }
     _objModel.makeSection(_VS);
     fireChange();
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Phi</em> parameter.
    *   This corresponds to Pitch or rotation around the X axis.
    *   @param rads the rotation angle in radians.
    */
   public void setPhi(double rads) {
      try {
        _obj.Wlz3DViewSetPhi(_VS, rads);
        _obj.WlzInit3DViewStruct(_VS,_obj);
      } catch (WlzException e) {
        System.err.println("setPhi");
        System.err.println(e);
      }
      _objModel.makeSection(_VS);
      fireChange();
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Zeta</em> parameter.
    *   This corresponds to Roll or rotation around the Y axis.
    *   @param rads the rotation angle in radians.
    */
   public void setZeta(double rads) {
      try {
        _obj.Wlz3DViewSetZeta(_VS, rads);
        _obj.WlzInit3DViewStruct(_VS,_obj);
      } catch (WlzException e) {
        System.err.println("setZeta");
        System.err.println(e);
      }
      _objModel.makeSection(_VS);
      fireChange();
   }

//----------------------------------------------
   // fixed line rotation
   /**
    *   Sets the ViewStructModel's <em>Psi</em> parameter.
    *   Note that this is <b>not</b> a View Structure parameter,
    *   it is only used by the ViewStructModel to calculate
    *   the change in Theta Phi and Zeta required to produce
    *   a rotation of Psi around the fixed line..
    *   @param rads the rotation angle in radians.
    */
   public void setPsi(double rads) {

      double phi = 0.0;
      double theta = 0.0;

      double phiDeg = 0.0;
      double thetaDeg = 0.0;
      double zetaDeg = 0.0;

      phi = Math.acos(Math.cos(rads)*_norm2[2] + Math.sin(rads)*_norm3[2]);
      if(Math.abs(Math.sin(phi)) < .01 ) {
         theta = 0.0;
      } else {
         theta = Math.atan2(Math.cos(rads)*_norm2[1] + Math.sin(rads)*_norm3[1],
                       Math.cos(rads)*_norm2[0] + Math.sin(rads)*_norm3[0]);
         if( theta < 0.0 ){
            theta += 2 * Math.PI;
         }
      }

      /*
      System.out.println("ViewStructModel.setPsi:");
      System.out.println("norm2 = "+Double.toString(_norm2[0])+", "+
                                    Double.toString(_norm2[1])+", "+
				    Double.toString(_norm2[2]));
      System.out.println("norm3 = "+Double.toString(_norm3[0])+", "+
                                    Double.toString(_norm3[1])+", "+
				    Double.toString(_norm3[2]));
      */
      phiDeg = phi * 180.0 / Math.PI;
      thetaDeg = theta * 180.0 / Math.PI;
      /*
      System.out.println("phi = "+ Double.toString(phiDeg));
      System.out.println("theta = "+ Double.toString(thetaDeg));
      */


      try {
        _obj.Wlz3DViewSetPhi(_VS, phi);
        _obj.Wlz3DViewSetTheta(_VS, theta);
        _obj.WlzInit3DViewStruct(_VS,_obj);
      } catch (WlzException e) {
        System.err.println("setPsi");
        System.err.println(e);
      }
      _objModel.makeSection(_VS);
      fireChange();
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Dist</em> parameter.
    *   @param dist orthogonal distance from the <em>fixed point</em>
    *   to the plane of the current section.
    */
   public void setDist(double dist) {

     try {
       _obj.Wlz3DViewSetDist(_VS, dist);
       _obj.WlzInit3DViewStruct(_VS,_obj);
     } catch (WlzException e) {
       System.err.println("setDist");
       System.err.println(e);
     }
     _objModel.makeSection(_VS);
     fireChange();
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>Fixed Point</em>.
    *   @param x x coordinate of the fixed point.
    *   @param y y coordinate of the fixed point.
    *   @param z z coordinate of the fixed point.
    */
   public void setFixedPoint(double x, double y,  double z) {
     try {
       _obj.Wlz3DViewSetFixed(_VS, x,y,z);
       _obj.Wlz3DViewSetDist(_VS, 0.0);
       _obj.WlzInit3DViewStruct(_VS,_obj);
     } catch (WlzException e) {
       System.err.println("setFixedPoint");
       System.err.println(e);
     }
     _objModel.makeSection(_VS);
     fireChange();
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>2nd Fixed Point</em>.
    *   This defines a <em>Fixed Line</em> for arbitrary rotation.
    *   @param x x coordinate of the 2nd fixed point.
    *   @param y y coordinate of the 2nd fixed point.
    *   @param z z coordinate of the 2nd fixed point.
    */
   public void setAxisPoint(double x, double y, double z) {

     try {
       _obj.Wlz3DViewSetFixed2(_VS, x,y,z);
       _obj.WlzInit3DViewStruct(_VS,_obj);
     } catch (WlzException e) {
       System.err.println("setAxisPoint");
       System.err.println(e);
     }
   }

//----------------------------------------------
   /**
    *   Sets the View Structure's <em>FixedLineAngle</em> parameter.
    *   This is used in <em>Fixed Line</em view mode.
    *   @param angle the angle (in radians) of the
    *   <em>fixed line</em> on a 2D section.
    */
   public void setFixedLineAngle(double angle) {
     try {
       _obj.Wlz3DViewSetFixedLineAngle(_VS, angle);
     } catch (Exception e) {
       System.err.println(e);
     }
   }

//----------------------------------------------
   /**
    *   Sets the ViewStructModel's <em>Norm2</em> parameter.
    *   Note that this is <b>not</b> a View Structure parameter,
    *   it is only used by the ViewStructModel to calculate
    *   the change in Theta Phi and Zeta required to produce
    *   a rotation of Psi around the fixed line..
    *   @param x x value of the norm2 vector.
    *   @param y y value of the norm2 vector.
    *   @param z z value of the norm2 vector.
    */
   public void setNorm2(double x, double y, double z) {
     _norm2[0] = x;
     _norm2[1] = y;
     _norm2[2] = z;
   }

//----------------------------------------------
   /**
    *   Sets the ViewStructModel's <em>Norm3</em> parameter.
    *   Note that this is <b>not</b> a View Structure parameter,
    *   it is only used by the ViewStructModel to calculate
    *   the change in Theta Phi and Zeta required to produce
    *   a rotation of Psi around the fixed line..
    *   @param x x value of the norm3 vector.
    *   @param y y value of the norm3 vector.
    *   @param z z value of the norm3 vector.
    */
   public void setNorm3(double x, double y, double z) {
     _norm3[0] = x;
     _norm3[1] = y;
     _norm3[2] = z;
   }

//==============================================
//==============================================
   /**
    *   Gets the View Structure's <em>Phi</em> parameter.
    *   @param phi a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    */
   public void getPhi(double[] phi) {
     try {
       _obj.Wlz3DViewGetPhi(_VS, phi);
     } catch (WlzException e) {
       System.err.println("getPhi");
       System.err.println(e);
     }
   }

   public double getPhi(){
     double[] phi = new double[1];
     phi[0] = 0.0;
     try {
       _obj.Wlz3DViewGetPhi(_VS, phi);
     } catch (WlzException e) {
       System.err.println("getPhi");
       System.err.println(e);
     }

     return phi[0];
   }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>Theta</em> parameter.
    *   @param theta a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    */
   public void getTheta(double[] theta) {
     try {
       _obj.Wlz3DViewGetTheta(_VS, theta);
     } catch (WlzException e) {
       System.err.println("getTheta");
       System.err.println(e);
     }
   }

   public double getTheta(){
     double theta[] = new double[1];
     theta[0] = 0.0;
     try {
        _obj.Wlz3DViewGetZeta(_VS, theta);
      } catch (WlzException e) {
        System.err.println("getZeta");
        System.err.println(e);
      }
     return theta[0];
   }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>Zeta</em> parameter.
    *   @param zeta a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    */
   public void getZeta(double[] zeta) {
      try {
        _obj.Wlz3DViewGetZeta(_VS, zeta);
      } catch (WlzException e) {
        System.err.println("getZeta");
        System.err.println(e);
      }
   }

   public double getZeta(){
     double zeta[] = new double[1];
     zeta[0] = 0.0;
     try {
        _obj.Wlz3DViewGetZeta(_VS, zeta);
      } catch (WlzException e) {
        System.err.println("getZeta");
        System.err.println(e);
      }
     return zeta[0];
   }
//----------------------------------------------
   /**
    *   Gets the View Structure's <em>Dist</em> parameter.
    *   @param dist a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    */
   public void getDist(double[] dist) {
      try {
        _obj.Wlz3DViewGetDist(_VS, dist);
      } catch (WlzException e) {
        System.err.println("getDist");
        System.err.println(e);
      }
   }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>Fixed Point</em>.
    *   @param x a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param y a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param z a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    */
   public void getFixedPoint(double[] x, double[] y, double[] z) {
      try {
        _obj.Wlz3DViewGetFixed(_VS, x,y,z);
      } catch (WlzException e) {
        System.err.println("getFixedPoint");
        System.err.println(e);
      }
   }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>Fixed Point</em>.
    *   @param x a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param y a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param z a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @return a 3 element array holding the
    *   <em>fixed point</em> coordinates.
    */
    public double[] getFixedPoint() {
      double[] x = new double[1];
      double[] y = new double[1];
      double[] z = new double[1];
      double[] fp = new double[3];

      try {
        _obj.Wlz3DViewGetFixed(_VS, x, y, z);
      }
      catch (WlzException e) {
        System.err.println("getFixedPoint");
        System.err.println(e);
      }

      fp[0] = x[0];
      fp[1] = y[0];
      fp[2] = z[0];

      return fp;
    }

//----------------------------------------------
   /**
    *   Gets the initial <em>fixed point</em> from
    *   the ViewStructModel.
    *   @return a 3 element array holding the initial
    *   <em>fixed point</em> coordinates.
    */
   public double[] getInitialFixedPoint() {
      return _fpInitial;
   }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>2nd Fixed Point</em>.
    *   @param x a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param y a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param z a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    */
   public void getAxisPoint(double[] x, double[] y, double[] z) {
      try {
        _obj.Wlz3DViewGetFixed2(_VS, x,y,z);
      } catch (WlzException e) {
        System.err.println("getAxisPoint");
        System.err.println(e);
      }
   }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>2nd Fixed Point</em>.
    *   @param x a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param y a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @param z a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    *   @return a 3 element array holding the
    *   2nd <em>fixed point</em> coordinates.
    */
    public double[] getAxisPoint() {
      double[] x = new double[1];
      double[] y = new double[1];
      double[] z = new double[1];
      double[] ap = new double[3];

      try {
        _obj.Wlz3DViewGetFixed2(_VS, x, y, z);
      }
      catch (WlzException e) {
        System.err.println("getAxisPoint");
        System.err.println(e);
      }

      ap[0] = x[0];
      ap[1] = y[0];
      ap[2] = z[0];

      return ap;
    }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>ViewMode</em> parameter.
    *   @return the description of the ViewMode returned by the View Structure.
    *   This will be one of
    *   <ul><li>WLZ_STATUE_MODE</li>
    *   <li>WLZ_UP_IS_UP_MODE  (corresponding to <em>up is up</em>)</li>
    *   <li>WLZ_FIXED_LINE_MODE  (corresponding to <em>fixed line</em>)</li>
    *   <li>WLZ_ZERO_ZETA_MODE</li>
    *   <li>WLZ_ZETA_MODE  (corresponding to <em>absolute</em>)</li></ul>
    */
   public String getViewMode() {
      String ret = "";
      iValArray[0] = 0;

      try {
        _obj.Wlz3DViewGetViewMode(_VS, iValArray);
      } catch (WlzException e) {
        System.err.println("getViewMode");
        System.err.println(e);
      }

      if(iValArray[0] == 0) {
        ret = "WLZ_STATUE_MODE\n";
      } else if(iValArray[0] == 1) {
        ret = "WLZ_UP_IS_UP_MODE\n";
      } else if(iValArray[0] == 2) {
        ret = "WLZ_FIXED_LINE_MODE\n";
      } else if(iValArray[0] == 3) {
        ret = "WLZ_ZERO_ZETA_MODE\n";
      } else if(iValArray[0] == 4) {
        ret = "WLZ_ZETA_MODE\n";
      } else {
        ret = "unrecognised view mode\n";
      }
      return ret;
   }

//----------------------------------------------
   /**
    *   Gets the View Structure's <em>FixedLineAngle</em> parameter.
    *   This is the angle (in radians) of the
    *   <em>fixed line</em> on a 2D section.
    *   @param arg a single element array used to retrieve
    *   the value of a <em>pointer to double</em> from
    *   a native C function.
    */
   public void getFixedLineAngle(double[] ang) {

     try {
       _obj.Wlz3DViewGetFixedLineAngle(_VS, ang);
     } catch (Exception e) {
       System.err.println(e);
     }
   }

//----------------------------------------------
   /**
    *   Gets the <em>inverse affine transform</em>
    *   corresponding to the current View Structure.
    *   @return the inverse affine transform.
    */
   public WlzAffineTransform getInverseTransform() {

     WlzAffineTransform ret = null;
     WlzAffineTransform trans = null;

     try {
       trans = _obj.WlzBndGetAffineTransform(_VS);
       if(trans == null) return null;
       ret = _obj.WlzAffineTransformInverse(trans);
     } catch (Exception e) {
       System.err.println(e);
     }

     return ret;
   }


//----------------------------------------------
   /**
    *   Gets the current View Structure.
    *   @return the current View Structure.
    */
   public WlzThreeDViewStruct getViewStruct() {
      return _VS;
   }

//----------------------------------------------
   /**
    *   Gets the <em>transformation matrix</em>
    *   corresponding to the given affine transform.
    *   @return a 3 dimensional array used to retrieve
    *   the value of a <em>pointer to transformation matrix</em>
    *   from a native C function.
    *   The first dimension of the array has a single element
    *   which is the transformation matrix.
    */
   public double[][][] getTransMatrix(WlzAffineTransform trans) {

     double mat[][][] = new double[1][][];
     mat[0] = null;

     WlzIVertex2 size[] = new WlzIVertex2[1];
     size[0] = null;

     try {
        _obj.WlzBndGetTransformMatrix(size, mat, trans);
     }
     catch(WlzException e) {
        System.out.println("WlzErr getting transform matrix");
        System.out.println(e.getMessage());
     }

     return mat;
   }

//----------------------------------------------
   /**
    *   Gets the ViewStructModel's <em>Norm2</em> vector.
    *   @return the _norm2 vector.
    */
   public double[] getNorm2() {
     return _norm2;
   }

//----------------------------------------------
   /**
    *   Gets the ViewStructModel's <em>Norm3</em> vector.
    *   @return the _norm3 vector.
    */
   public double[] getNorm3() {
     return _norm3;
   }

//----------------------------------------------
   /**
    *   Utility to print double array's for debugging..
    *   @param arr the array to print.
    */
   public void printArray(double[] arr) {
      int len = arr.length;

      for(int i=0; i<len; i++) {
         System.out.println("element "+i+" = "+arr[i]);
      }

   }

//----------------------------------------------
   /**
    *   Utility to print the current View Structure for debugging..
    */
   public void printViewStructure() {
      System.out.println("VIEW STRUCTURE");
      try {
        _obj.Wlz3DViewGetTheta(_VS, dValArray);
        System.out.println("theta = " + dValArray[0]);
        dValArray[0] = 0.0d;
        _obj.Wlz3DViewGetPhi(_VS, dValArray);
        System.out.println("phi = " + dValArray[0]);
        dValArray[0] = 0.0d;
        _obj.Wlz3DViewGetZeta(_VS, dValArray);
        System.out.println("zeta = " + dValArray[0]);
        dValArray[0] = 0.0d;
        _obj.Wlz3DViewGetDist(_VS, dValArray);
        System.out.println("dist = " + dValArray[0]);
        _obj.Wlz3DViewGetFixed(_VS, dValArray1, dValArray2, dValArray3);
        System.out.println("fixed = " + dValArray1[0] + ", " +
                                        dValArray2[0] + ", " + dValArray3[0]);
        dValArray[0] = 0.0d;
        _obj.Wlz3DViewGetScale(_VS, dValArray);
        System.out.println("scale = " + dValArray[0]);
        iValArray[0] = 0;
        _obj.Wlz3DViewGetViewMode(_VS, iValArray);
        if(iValArray[0] == 0) {
          System.out.println("view mode = WLZ_STATUE_MODE\n");
        } else if(iValArray[0] == 1) {
          System.out.println("view mode = WLZ_UP_IS_UP_MODE\n");
        } else if(iValArray[0] == 2) {
          System.out.println("view mode = WLZ_FIXED_LINE_MODE\n");
        } else if(iValArray[0] == 3) {
          System.out.println("view mode = WLZ_ZERO_ZETA_MODE\n");
        } else if(iValArray[0] == 4) {
          System.out.println("view mode = WLZ_ZETA_MODE\n");
        } else {
          System.out.println("unrecognised view mode\n");
        }
      }
      catch (WlzException e) {
        System.err.println("printViewStructure");
        System.err.println(e);
      }
      System.out.println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   } // printViewStructure()

//-------------------------------------------------------------
// handle all objects that are interested in changes
//-------------------------------------------------------------

   // keep track of all the listeners to this model
   /**
    *   A list of ChangeListeners which are
    *   listening for events fired from the ViewStructModel.
    */
   protected EventListenerList changeListeners = new EventListenerList();

   // add a listener to the register
   /**
    *   Registers a ChangeListener
    *   with the EventListenerList.
    *   @param x an Event handler implementing ChangeListener
    */
   public void addChangeListener(ChangeListener x) {
      changeListeners.add (ChangeListener.class, x);
   }

   // remove a listener from the register
   /**
    *   Removes a ChangeListener
    *   from the EventListenerList.
    *   @param x an Event handler implementing ChangeListener
    */
   public void removeChangeListener(ChangeListener x) {
      changeListeners.remove (ChangeListener.class, x);
   }

   /**  An event to be fired from the ViewStructModel */
   private ChangeEvent ce;

   /**  A local copy of the list of Change Listeners */
   private Object[] listeners;

   /**  One of the list of Change Listeners */
   private ChangeListener cl;

   /**  A counter */
   private int ijk = 0;

   /**
    *   Fires a ChangeEvent from the ViewStructModel.
    */
   public void fireChange() {
   //   if(_debug) System.out.println("VSModel.fireChange()");
      // Create the event:
      ce = new ChangeEvent(this);
      // Get the listener list
      listeners = changeListeners.getListenerList();
      // Process the listeners last to first
      // List is in pairs, Class and instance
      for (ijk = listeners.length-2; ijk >= 0; ijk -= 2) {
        if (listeners[ijk] == ChangeListener.class) {
          cl = (ChangeListener)listeners[ijk + 1];
          cl.stateChanged(ce);
        }
      }
   } // fireChange

//-------------------------------------------------------------
} // class ViewStruct
