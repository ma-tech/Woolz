/************************************************************************
* Project:      Java Woolz
* Title:        WlzFBox3.java
* Date:         January 1999
* Purpose:      Java object to mirror the Woolz WlzFBox2 structure.
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

public class WlzFBox3 extends WlzFBox2 implements Cloneable
{
  public float	zMin;
  public float	zMax;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzFBox3()
  {
    zMin = 0.0;
    zMax = 0.0;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     xmin		minimum box x value
  * @param:     ymin		minimum box y value
  * @param:     zmin		minimum box z value
  * @param:     xmax		maximum box x value
  * @param:     ymax		maximum box y value
  * @param:     zmax		maximum box z value
  **********************************************************************/
  public	WlzFBox3(float xmin, float ymin, float zmin,
  			 float xmax, float ymax, float zmax)
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
    return(new WlzFBox3(xMin, yMin, zMin, xMax, yMax, zMax));
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
      isEqual = (Math.abs(xMin -
      			  ((WlzFBox3 )other).xMin) < Double.MIN_VALUE) &&
                (Math.abs(yMin -
			  ((WlzFBox3 )other).yMin) < Double.MIN_VALUE) &&
                (Math.abs(zMin -
			  ((WlzFBox3 )other).zMin) < Double.MIN_VALUE) &&
                (Math.abs(xMax -
			  ((WlzFBox3 )other).xMax) < Double.MIN_VALUE) &&
                (Math.abs(yMax -
			  ((WlzFBox3 )other).yMax) < Double.MIN_VALUE) &&
                (Math.abs(zMax -
			  ((WlzFBox3 )other).zMax) < Double.MIN_VALUE);
    }
    return(isEqual);
  }
}
