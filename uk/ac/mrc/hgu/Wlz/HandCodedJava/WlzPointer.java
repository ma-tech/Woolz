/************************************************************************
* Project:      Java Woolz
* Title:        WlzPointer.java
* Date:         January 1999
* Purpose:      Java binding for Woolz native pointers.
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

public class WlzPointer extends WlzNative
{
  /**********************************************************************
  * Purpose:    Indicates whether some other object is "equal to" this
  *		Woolz pointer.
  * @return     true if this object is the same as the given object,
  *		otherwise false.
  * @param:     obj		the given object for comparison.
  **********************************************************************/
  public boolean equals(Object other)
  {
    boolean	isEqual;

    if(other == null)
    {
      isEqual = false;
    }
    else if(other == this)
    {
      isEqual = true;
    }
    else if(WlzPointer.class != other.getClass())
    {
      isEqual = false;
    }
    else
    {
      isEqual = value == ((WlzPointer )other).value;
    }
    return(isEqual);
  }

  /**********************************************************************
  * Purpose:    Computes a hashcode for this pointer.
  * @return     a hash code value for the Woolz pointer.
  * @param:     void
  **********************************************************************/
  public int 	hashCode()
  {
    return((int )((value >>> 32) + (value & 0xFFFFFFFF)));
  }

  /**********************************************************************
  * Purpose:    Tests the pointer to see if it is null.
  * @return     true if the pointer is null.
  * @param:     void
  **********************************************************************/
  public boolean isNull()
  {
    return(value == 0);
  }

  /**********************************************************************
  * Purpose:    Makes the pointer a null pointer.
  * @return     void
  * @param:     void
  **********************************************************************/
  public void	clear()
  {
    value = 0;
  }
}
