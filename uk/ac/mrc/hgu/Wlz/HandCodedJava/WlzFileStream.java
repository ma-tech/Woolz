/************************************************************************
* Project:      Java Woolz
* Title:        WlzFileStream.java
* Date:         January 1999
* Purpose:      Java binding for Woolz native file input.
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
import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzFileStream extends WlzPointer
{
  private native long	 JWlzOpen(String name, String rw) throws IOException;
  private native void	 JWlzClose(long value) throws IOException;

  /**********************************************************************
  * Purpose:    Closes the native file input stream.
  * @return     void
  * @param:     void
  **********************************************************************/
  protected void finalize()
  throws IOException
  {
    close();
  }

  /**********************************************************************
  * Purpose:    Creates a new native file input stream.
  * @param:     name		the system-dependent file name.
  **********************************************************************/
  public WlzFileStream(String name, String rw)
  throws IOException
  {
    value = JWlzOpen(name, rw);
  }

  /**********************************************************************
  * Purpose:    Closes the native file input stream.
  * @return     void
  * @param:     void
  **********************************************************************/
  public void	close()
  throws IOException
  {
    JWlzClose(value);
  }
}
