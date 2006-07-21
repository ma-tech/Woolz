/************************************************************************
* Project:      Java Woolz
* Title:        WlzFBox2.java
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

public class WlzFBox2 extends WlzBase implements Cloneable
{
  public float	xMin;
  public float	xMax;
  public float	yMin;
  public float	yMax;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzFBox2()
  {
    xMin = 0.0f;
    xMax = 0.0f;
    yMin = 0.0f;
    yMax = 0.0f;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     xmin		minimum box x value
  * @param:     ymin		minimum box y value
  * @param:     xmax		maximum box x value
  * @param:     ymax		maximum box y value
  **********************************************************************/
  public	WlzFBox2(float xmin, float ymin, float xmax, float ymax)
  {
    xMin = xmin;
    xMax = xmax;
    yMin = ymin;
    yMax = ymax;
  }

  /**********************************************************************
  * Purpose:	Constructs a new Java (2D) rectangle which has the same
  *		coordinates as this Java Woolz box.
  * @param:	void
  **********************************************************************/
  public Rectangle2D toRectangle()
  {
    return(new Rectangle2D.Double(xMin, yMin, xMax - xMin, yMax - yMin));
  }

  /**********************************************************************
  * Purpose:	Implements cloning
  * @return:	Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object	clone()
  {
    return(new WlzFBox2(xMin, yMin, xMax, yMax));
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
      			  ((WlzFBox2 )other).xMin) < Double.MIN_VALUE) &&
                (Math.abs(yMin -
			  ((WlzFBox2 )other).yMin) < Double.MIN_VALUE) &&
                (Math.abs(xMax -
			  ((WlzFBox2 )other).xMax) < Double.MIN_VALUE) &&
                (Math.abs(yMax -
			  ((WlzFBox2 )other).yMax) < Double.MIN_VALUE);
    }
    return(isEqual);
  }
}
