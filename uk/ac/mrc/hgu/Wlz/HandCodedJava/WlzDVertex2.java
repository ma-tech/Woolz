/************************************************************************
* Project:      Java Woolz
* Title:        WlzDVertex2.java
* Date:         January 1999
* Purpose:      Java object to mirror the Woolz WlzDVertex2 structure.
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

public class WlzDVertex2 extends WlzBase implements Cloneable
{
  public double	vtX;
  public double	vtY;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzDVertex2()
  {
    vtX = 0.0;
    vtY = 0.0;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     x		the x coordinate
  * @param:     y		the y coordinate
  **********************************************************************/
  public	WlzDVertex2(double x, double y)
  {
    vtX = x;
    vtY = y;
  }

  /**********************************************************************
  * Purpose:	Constructs a new Java (2D) point which has the same
  *		coordinates as this Java Woolz vertex.
  * @return	new Java (2D) point
  * @param:	void
  **********************************************************************/
  public Point2D toPoint()
  {
    return(new Point2D.Double(vtX, vtY));
  }

  /**********************************************************************
  * Purpose:    Implements cloning.
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzDVertex2(vtX, vtY));
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
			  ((WlzDVertex2 )other).vtY) < Double.MIN_VALUE) &&
                (Math.abs(vtY -
			  ((WlzDVertex2 )other).vtY) < Double.MIN_VALUE);
    }
    return(isEqual);
  }
}
