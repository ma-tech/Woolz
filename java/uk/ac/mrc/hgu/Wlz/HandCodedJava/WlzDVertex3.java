/************************************************************************
* Project:      Java Woolz
* Title:        WlzDVertex3.java
* Date:         January 1999
* Purpose:      Java object to mirror the Woolz WlzDVertex3 structure.
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

import java.lang.*;
import java.awt.geom.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzDVertex3 extends WlzDVertex2 implements Cloneable
{
  public double	vtZ;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzDVertex3()
  {
    vtZ = 0.0;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     x		the x coordinate
  * @param:     y		the y coordinate
  * @param:     z		the z coordinate
  **********************************************************************/
  public	WlzDVertex3(double x, double y, double z)
  {
    vtX = x;
    vtY = y;
    vtZ = z;
  }

  /**********************************************************************
  * Purpose:    Implements cloning * @param:     void
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzDVertex3(vtX, vtY, vtZ));
  }

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
      isEqual = (Math.abs(vtX -
			  ((WlzDVertex3 )other).vtY) < Double.MIN_VALUE) &&
                (Math.abs(vtY -
			  ((WlzDVertex3 )other).vtY) < Double.MIN_VALUE) &&
                (Math.abs(vtZ -
			  ((WlzDVertex3 )other).vtZ) < Double.MIN_VALUE);
    }
    return(isEqual);
  }
}
