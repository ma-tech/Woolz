/************************************************************************
* Project:      Java Woolz
* Title:        WlzException.java
* Date:         January 1999
* Purpose:      Java exception for Woolz errors.
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

public class WlzException extends Exception
{
  /**********************************************************************
  * Purpose:    Constructor, just to make it explicit.
  * @param:     void
  **********************************************************************/
  public WlzException()
  {
    super();
  }

  /**********************************************************************
  * Purpose:    Constructor, with detail message.
  * @param:     message		the detail message.
  **********************************************************************/
  public WlzException(String message)
  {
    super(message);
  }
}
