/************************************************************************
* Project:      Java Woolz
* Title:        WlzPixelP.java
* Date:         January 1999
* Purpose:      Java binding for Woolz grey pointer structure.
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

public class WlzPixelP extends WlzGreyP
{
  protected int	type;

  /**********************************************************************
  * Purpose:    Get the grey value type
  * @return     int		the grey value type of the pixel
  * @param:     void
  **********************************************************************/
  public int	getType()
  {
    return(type);
  }
}
