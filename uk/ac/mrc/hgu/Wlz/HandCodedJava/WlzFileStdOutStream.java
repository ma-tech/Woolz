/************************************************************************
* Project:      Java Woolz
* Title:        WlzFileStdOutStream.java
* Date:         January 1999
* Purpose:      Java binding for Woolz native file output.
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

public class WlzFileStdOutStream extends WlzFileStream
{
  /**********************************************************************
  * Purpose:    Creates a new native standard output file stream.
  * @param:     void
  **********************************************************************/
  public WlzFileStdOutStream()
  throws IOException
  {
    super("stdout", "");
  }
}
