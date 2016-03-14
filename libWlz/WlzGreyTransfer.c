#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyTransfer_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreyTransfer.c
* \author       Richard Baldock, Bill Hill
* \date         September 2002
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
* \brief	Transfers grey values from a source to a destination
* 		object. The source object grey values are set in the
* 		intersection domain between source and destination.
* 		Destination domain and the destination values outside
* 		of the intersection are unchanged.
* \ingroup	WlzValuesUtils
*/

#include <stdio.h>
#include <string.h>
#include <Wlz.h>

static WlzErrorNum		WlzGreyTransfer2D(
				  WlzObject *dObj,
				  WlzObject *sObj);


/*!
* \return	New object or NULL on error.
* \ingroup	WlzValuesUtils
* \brief 	Transfers grey values from the source object to the
*               destination object within the intersection of the source
*               and destination. Grey values within the destination
*               object outside of the source object are unchanged.
*               It is an error if either object has a different dimension
*               or grey value type, except for when either is an empty
*               object.
* \param	dObj			Destination object which may be
* 					empty, but otherwise should be of the
* 					same dimension as the source object
* 					with valid values..
* \param	sObj			Source object which if not empty must
* 					have both a valid domain and valid
* 					values.
* \param	inplace			Overwrite the destination object's
* 					values if non zero.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzGreyTransfer(
				  WlzObject *dObj,
				  WlzObject *sObj,
				  int inplace,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(WlzIsEmpty(dObj, NULL))
  {
    rObj = WlzMakeEmpty(&errNum);
  }
  else if(WlzIsEmpty(sObj, NULL))
  {
    rObj = WlzMakeMain(dObj->type, dObj->domain, dObj->values,
                       dObj->plist, NULL, &errNum);
  }
  else if(dObj->type != sObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((dObj->domain.core == NULL) || (sObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dObj->values.core == NULL) || (sObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(sObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ: /* FALLTHROUGH */
        {
	  WlzObject	*rIObj = NULL;

	  rIObj = WlzIntersect2(dObj, sObj, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    rObj = (inplace)?
		   WlzMakeMain(dObj->type, dObj->domain, dObj->values,
			       dObj->plist, NULL, &errNum):
		   WlzCopyObject(dObj, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(sObj->type == WLZ_2D_DOMAINOBJ)
	    {
	      WlzObject *sIObj;

	      rIObj->values = WlzAssignValues(rObj->values, NULL);
	      sIObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	                          rIObj->domain, sObj->values,
				  NULL, NULL, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        errNum = WlzGreyTransfer2D(rIObj, sIObj);
	      }
	      (void )WlzFreeObj(sIObj);
	    }
	    else /* sObj->type == WLZ_3D_DOMAINOBJ */
	    {
	      int	p,
			rTiled,
			sTiled,
	      		nPlanes;

	      rTiled = WlzGreyTableIsTiled(rObj->values.core->type);
	      sTiled = WlzGreyTableIsTiled(sObj->values.core->type);
	      nPlanes = rIObj->domain.p->lastpl - rIObj->domain.p->plane1 + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	      for(p = 0; p < nPlanes; ++p)
	      {
		if(errNum == WLZ_ERR_NONE)
		{
		  int	     pln;
		  WlzDomain  dom;
		  WlzValues  val;
		  WlzObject  *rIObj2D = NULL,
			     *sIObj2D = NULL;
                  WlzErrorNum errNum2D = WLZ_ERR_NONE;

		  pln = p + rIObj->domain.p->plane1;
		  dom = rIObj->domain.p->domains[p];
		  val = (rTiled)?
		        rObj->values:
		        rObj->values.vox->values[pln -
			                         rObj->values.vox->plane1];
		  rIObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ,
		                        dom, val, NULL, NULL, &errNum2D);
		  if(errNum2D == WLZ_ERR_NONE)
		  {
		    val = (sTiled)?
		          sObj->values:
			  sObj->values.vox->values[pln -
			                           sObj->values.vox->plane1];
		    sIObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ,
		                          dom, val, NULL, NULL, &errNum2D);
		  }
		  if(errNum2D == WLZ_ERR_NONE)
		  {
		    errNum2D = WlzGreyTransfer2D(rIObj2D, sIObj2D);
		  }
		  (void )WlzFreeObj(rIObj2D);
		  (void )WlzFreeObj(sIObj2D);
#ifdef _OPENMP
#pragma omp critical
		  {
#endif
		    if((errNum == WLZ_ERR_NONE) && (errNum2D != WLZ_ERR_NONE))
		    {
		      errNum = errNum2D;
		    }
#ifdef _OPENMP
		  }
#endif
		}
	      }
	    }
	  }
	  (void )WlzFreeObj(rIObj);
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief 	Transfers grey values from the 2D source object to the
*               2D destination object where both share the same domain.
*               The objects are assumed valid 2D domain objects from
*               WlzGreyTransfer().
* \param	dObj			Destination object.
* \param	sObj			Source object.
*/
static WlzErrorNum		WlzGreyTransfer2D(
				  WlzObject *dObj,
				  WlzObject *sObj)
{
  WlzGreyWSpace sGWSp,
  		dGWSp;
  WlzIntervalWSpace sIWSp,
  		    dIWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(sObj, &sIWSp, &sGWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(dObj, &dIWSp, &dGWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		bub;

    bub = (sGWSp.pixeltype == dGWSp.pixeltype) &&
          (sGWSp.pixeltype == WLZ_GREY_UBYTE);
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&sIWSp)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&dIWSp)) == WLZ_ERR_NONE))
    {
      if(bub)
      {
	/* Avoid switches and function calls for the most common case. */
	(void )memcpy(dGWSp.u_grintptr.ubp, sGWSp.u_grintptr.ubp,
		      sIWSp.colrmn * sizeof(WlzUByte));
      }
      else
      {
	WlzValueCopyGreyToGrey(dGWSp.u_grintptr, 0, dGWSp.pixeltype,
			       sGWSp.u_grintptr, 0, sGWSp.pixeltype,
			       sIWSp.colrmn);
      }
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    (void )WlzEndGreyScan(&sIWSp, &sGWSp);
    (void )WlzEndGreyScan(&dIWSp, &dGWSp);
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}
