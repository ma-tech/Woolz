#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRGBAGreyStats_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzRGBAGreyStats.c
* \author       Richard Baldock
* \date         June 2005
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
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
* \brief	Calculates simple quick statistics for a domain object
* 		with RGBA values.
* \ingroup	WlzFeatures
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return       Object area or -1 on error.
* \ingroup      WlzFeatures
* \brief        Calculates simple quick statistics for the domain
*               object with RGBA values. Each component has its
*               statictics computed and entered into the four
*               double[4] arrays.
* \param        srcObj                  Given object.
* \param	colSpc			Colour space.
* \param        dstGType                Pointer for grey type.
* \param        dstMin                  Array for the 4 minimum value.
* \param        dstMax                  Array for the 4 maximum value.
* \param        dstSum                  Array for the 4 sum of values.
* \param        dstSumSq                Array for the 4 sum of squares of
*                                       values.
* \param	dstMean			Array for the 4 mean values.
* \param	dstStdDev		Array for the 4 standard deviation
* 					values.
* \param        dstErr                  Destination pointer for error, may
*                                       be NULL.
*/
int WlzRGBAGreyStats(
  WlzObject	*srcObj,
  WlzRGBAColorSpace	colSpc,
  WlzGreyType	*dstGType,
  double 	*dstMin,
  double 	*dstMax,
  double 	*dstSum,
  double 	*dstSumSq,
  double 	*dstMean,
  double 	*dstStdDev,
  WlzErrorNum 	*dstErr)
{
  int		area = 0;
  WlzCompoundArray	*cmpnd;
  int		i;
  WlzGreyType	gType;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* all object checks done bt RGBAToCOmpound including
     grey-level type */
  if((cmpnd = WlzRGBAToCompound(srcObj, colSpc, &errNum)) != NULL){
    *dstGType = WLZ_GREY_RGBA;
    for(i=0; (i < 4) && (errNum == WLZ_ERR_NONE) ; i++){
      area = WlzGreyStats(cmpnd->o[i], &gType,
			  dstMin + i, dstMax + i,
			  dstSum + i, dstSumSq + i,
			  dstMean + i, dstStdDev + i,
			  &errNum);
    }
    WlzFreeObj((WlzObject *) cmpnd);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(area);
}
