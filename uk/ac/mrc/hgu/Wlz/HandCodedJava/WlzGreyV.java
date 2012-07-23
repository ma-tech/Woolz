/*!
* @file         WlzGreyV.java
* @author       Bill Hill
* @date         January 1999
* @version      $Id$
* @par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* @par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* @brief        Java binding for Woolz grey value structure.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;
import java.lang.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzGreyV extends WlzNative implements Cloneable
{
  // Version string.
  public static String ident = "Id$$";

  /*!
  * @brief      Constructor for a long grey value.
  * @param      v			The grey value.
  */
  public	WlzGreyV(long v)
  {
    setLongValue(v);
  }

  /*!
  * @brief      Constructor for a long grey value.
  * @param      v			The grey value.
  */
  public	WlzGreyV(int v)
  {
    setIntValue(v);
  }

  /*!
  * @brief      Constructor for a long grey value.
  * @param      v			The grey value.
  */
  public	WlzGreyV(short v)
  {
    setShortValue(v);
  }

  /*!
  * @brief      Constructor for a long grey value.
  * @param      v			The grey value.
  */
  public	WlzGreyV(byte v)
  {
    setByteValue(v);
  }

  /*!
  * @brief      Constructor for a long grey value.
  * @param      v			The grey value.
  */
  public	WlzGreyV(float v)
  {
    setFloatValue(v);
  }

  /*!
  * @return     Clone of this object.
  * @brief      Implements cloning.
  */
  public Object clone()
  {
    return(new WlzGreyV(value));
  }
}
