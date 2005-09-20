#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzMakeCompound.c
* \author       Richard Baldock
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Makes compound objects.
* \ingroup	WlzAllocation
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

/*!
* \return	New Woolz compound array object.
* \ingroup	WlzAllocation
* \brief	Makes a new Woolz compound array object.
* \param	type			Type of compound array which should
*					be either WLZ_COMPOUND_ARR_1 or
*					WLZ_COMPOUND_ARR_2.
* \param	mode			Action to be performed, which may have
*					the values:
*					<ul>
*					<li>
*					1: Allocate an empty array with
*					space for n objects.
*					</li>
*					<li>
* 					2: Array members are given in the
*					input object list, the list linked.
*					</li>
*					<li>
* 					3: Array members are given in the
*					input object list, objects assigned.
*					</li>
*					</ul>
* \param	n			Number of objects.
* \param	ol			Input object list.
* \param	otype			Woolz object type if objects are
*					type checked.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzCompoundArray *
WlzMakeCompoundArray(WlzObjectType	type,
		     int 		mode,
		     int 		n,
		     WlzObject 		**ol,
		     WlzObjectType	otype,
		     WlzErrorNum	*dstErr)
{
  WlzCompoundArray	*co=NULL;
  int			i;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  
  if (type != WLZ_COMPOUND_ARR_1 && type != WLZ_COMPOUND_ARR_2) {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  /*
   * If appropriate, check objects have suitable type
   */
  if ((errNum == WLZ_ERR_NONE) && (type == WLZ_COMPOUND_ARR_1)) {
    switch(mode) {

    case 1:
      break;

    case 2:
    case 3:
      for (i=0; i<n; i++) {
	if (ol[i] != NULL && ol[i]->type != otype) {
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
	}
      }
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* Allocate space and fill if appropriate
     Note no default case because the mode has been checked above */
  if( errNum == WLZ_ERR_NONE ){
    switch (mode) {

    case 1:
      if( (co = (WlzCompoundArray *)
	   AlcCalloc(1,sizeof(WlzCompoundArray)
		     +n*sizeof(WlzObject *))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      co->o = (WlzObject **) (co+1);
      break;

    case 2:
      if( (co = (WlzCompoundArray *)
	   AlcCalloc(1,sizeof(WlzCompoundArray))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      co->o = ol;
      break;

    case 3:
      if( (co = (WlzCompoundArray *)
	   AlcCalloc(1,sizeof(WlzCompoundArray)
		     +n*sizeof(WlzObject *))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      co->o = (WlzObject **) (co+1);
      for (i=0; i<n; i++) {
	co->o[i] = WlzAssignObject(ol[i], NULL);
      }
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    co->n = n;
    co->type = type;
    if (type == WLZ_COMPOUND_ARR_1){
      co->otype = otype;
    }
    co->linkcount = 0;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(co);
}
