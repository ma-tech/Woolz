/************************************************************************
* Project:      Java Woolz
* Title:        WlzNative.java
* Date:         January 1999
* Purpose:      Java binding for Woolz native types; ie: int, short,
*		byte, float, double and pointer.
* Copyright:	1997 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Maintenance:	Log changes below, with most recent at top of list.
* @author       Bill Hill (bill@hgu.mrc.ac.uk)
* @version 	MRC HGU %I%, %G%
************************************************************************/
package uk.ac.mrc.hgu.Wlz;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzNative extends WlzBase
{
  // Possible types held in the native value union.
  public static final int LONG = 0;
  public static final int INT = 1;
  public static final int SHORT = 2;
  public static final int BYTE = 3;
  public static final int FLOAT = 4;
  public static final int DOUBLE = 5;
  public static final int POINTER = 6;

  // The native value union.
  protected long value;

  // Native methods.
  private native long JWlzNativeUnionFromLong(long val);
  private native long JWlzNativeUnionFromInt(int val);
  private native long JWlzNativeUnionFromShort(short val);
  private native long JWlzNativeUnionFromByte(byte val);
  private native long JWlzNativeUnionFromFloat(float val);
  private native long JWlzNativeUnionFromDouble(double val);
  private native long JWlzNativeUnionToLong(long value);
  private native int JWlzNativeUnionToInt(long value);
  private native short JWlzNativeUnionToShort(long value);
  private native byte JWlzNativeUnionToByte(long value);
  private native float JWlzNativeUnionToFloat(long value);
  private native double JWlzNativeUnionToDouble(long value);

  /**********************************************************************
  * Purpose:    Constructor, just to make it explicit.
  * @param:     void
  **********************************************************************/
  public WlzNative()
  {
    value = 0;
  }

  /**********************************************************************
  * Purpose:    Constructor with initial value which is only seen from
  *		the JNI C side.
  * @param:     val			Initial value.
  **********************************************************************/
  private WlzNative(long v)
  {
    value = v;
  }

  /**********************************************************************
  * Purpose:    Set the value to the given long value.
  * @param:     val			Given value.
  **********************************************************************/
  protected void setLongValue(long val)
  {
    value = val;
    // value = JWlzNativeUnionFromLong(val);
  }

  /**********************************************************************
  * Purpose:    Set the value to the given int value.
  * @param:     val			Given value.
  **********************************************************************/
  protected void setIntValue(int val)
  {
    value = JWlzNativeUnionFromInt(val);
  }

  /**********************************************************************
  * Purpose:    Set the value to the given short value.
  * @param:     val			Given value.
  **********************************************************************/
  protected void setShortValue(short val)
  {
    value = JWlzNativeUnionFromShort(val);
  }

  /**********************************************************************
  * Purpose:    Set the value to the given byte value.
  * @param:     val			Given value.
  **********************************************************************/
  protected void setByteValue(byte val)
  {
    value = JWlzNativeUnionFromByte(val);
  }

  /**********************************************************************
  * Purpose:    Set the value to the given float value.
  * @param:     val			Given value.
  **********************************************************************/
  protected void setFloatValue(float val)
  {
    value = JWlzNativeUnionFromFloat(val);
  }

  /**********************************************************************
  * Purpose:    Set the value to the given double value.
  * @param:     val			Given value.
  **********************************************************************/
  protected void setDoubleValue(double val)
  {
    value = JWlzNativeUnionFromDouble(val);
  }

  /**********************************************************************
  * Purpose:    Get the long value.
  * @param:     void
  **********************************************************************/
  protected long getLongValue()
  {
    // return(JWlzNativeUnionToLong(value));
    return(value);
  }

  /**********************************************************************
  * Purpose:    Get the int value.
  * @param:     void
  **********************************************************************/
  protected int getIntValue()
  {
    return(JWlzNativeUnionToInt(value));
  }

  /**********************************************************************
  * Purpose:    Get the short value.
  * @param:     void
  **********************************************************************/
  protected short getShortValue()
  {
    return(JWlzNativeUnionToShort(value));
  }

  /**********************************************************************
  * Purpose:    Get the byte value.
  * @param:     void
  **********************************************************************/
  protected byte getByteValue()
  {
    return(JWlzNativeUnionToByte(value));
  }

  /**********************************************************************
  * Purpose:    Get the float value.
  * @param:     void
  **********************************************************************/
  protected float getFloatValue()
  {
    return(JWlzNativeUnionToFloat(value));
  }

  /**********************************************************************
  * Purpose:    Get the double value.
  * @param:     void
  **********************************************************************/
  protected double getDoubleValue()
  {
    return(JWlzNativeUnionToDouble(value));
  }

}
