/************************************************************************
* Project:      Java Woolz
* Title:        WlzIVertex2.java
* Date:         January 1999
* Purpose:      Java object to mirror the Woolz WlzIVertex2 structure.
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
import java.awt.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzIVertex2 extends WlzBase implements Cloneable
{
  public int	vtX;
  public int	vtY;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzIVertex2()
  {
    vtX = 0;
    vtY = 0;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     x		the x coordinate
  * @param:     y		the y coordinate
  **********************************************************************/
  public	WlzIVertex2(int x, int y)
  {
    vtX = x;
    vtY = y;
  }

  /**********************************************************************
  * Purpose:	Constructs a new Java (2D) point which has the same
  *		coordinates as this Java Woolz vertex.
  * @param:	void
  **********************************************************************/
  public Point toPoint()
  {
    return(new Point(vtX, vtY));
  }

  /**********************************************************************
  * Purpose:    Implements cloning.
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzIVertex2(vtX, vtY));
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
      isEqual = (vtX == ((WlzIVertex2 )other).vtX) &&
                (vtY == ((WlzIVertex2 )other).vtY);
    }
    return(isEqual);
  }
}
