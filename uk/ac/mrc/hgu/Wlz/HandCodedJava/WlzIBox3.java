/************************************************************************
* Project:      Java Woolz
* Title:        WlzIBox3.java
* Date:         January 1999
* Purpose:      Java object to mirror the Woolz WlzIBox3 structure.
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

public class WlzIBox3 extends WlzIBox2 implements Cloneable
{
  public int	zMin;
  public int	zMax;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzIBox3()
  {
    zMin = 0;
    zMax = 0;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     xmin		minimum box x value
  * @param:     ymin		minimum box y value
  * @param:     xmax		maximum box x value
  * @param:     ymax		maximum box y value
  **********************************************************************/
  public	WlzIBox3(int xmin, int ymin, int zmin,
  			 int xmax, int ymax, int zmax)
  {
    xMin = xmin;
    xMax = xmax;
    yMin = ymin;
    yMax = ymax;
    zMin = zmin;
    zMax = zmax;
  }

  /**********************************************************************
  * Purpose:    Implements cloning.
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzIBox3(xMin, yMin, zMin, xMax, yMax, zMax));
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
      isEqual = (xMin == ((WlzIBox3 )other).xMin) &&
                (yMin == ((WlzIBox3 )other).yMin) &&
                (zMin == ((WlzIBox3 )other).zMin) &&
                (xMax == ((WlzIBox3 )other).xMax) &&
                (yMax == ((WlzIBox3 )other).yMax) &&
                (zMax == ((WlzIBox3 )other).zMax);
    }
    return(isEqual);
  }
}
