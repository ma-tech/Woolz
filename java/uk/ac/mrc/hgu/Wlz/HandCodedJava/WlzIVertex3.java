/************************************************************************
* Project:      Java Woolz
* Title:        WlzIVertex3.java
* Date:         January 1999
* Purpose:      Java object to mirror the Woolz WlzIVertex3 structure.
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

public class WlzIVertex3 extends WlzIVertex2 implements Cloneable
{
  public int	vtZ;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzIVertex3()
  {
    vtZ = 0;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     x		the x coordinate
  * @param:     y		the y coordinate
  * @param:     z		the z coordinate
  **********************************************************************/
  public	WlzIVertex3(int x, int y, int z)
  {
    vtX = x;
    vtY = y;
    vtZ = z;
  }

  /**********************************************************************
  * Purpose:    Implements cloning.
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzIVertex3(vtX, vtY, vtZ));
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
      isEqual = (vtX == ((WlzIVertex3 )other).vtX) &&
                (vtY == ((WlzIVertex3 )other).vtY) &&
                (vtZ == ((WlzIVertex3 )other).vtZ);
    }
    return(isEqual);
  }
}
