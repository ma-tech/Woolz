#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzClipObjToBox_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzClipObjToBox.c
* \author       Bill Hill
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
* \brief	Functions for clipping the domain of either 2D or 3D
* 		domain objects so that they lie within the given axis aligned
* 		biunding box.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	New object with clipped domain or NULL on error.
* \brief	Clips the given object's domain so that it does not lie outside
* 		the given 2D axis aliagned clip box. Only 2D domain objects and
* 		empty objects may be clipped by this function.
* \param	srcObj			Given source object.
* \param	clipBox			Clip box.
* \param	dstErrNum		Destination pointer for error, may be
* 					NULL.
*/
WlzObject	*WlzClipObjToBox2D(WlzObject *srcObj, WlzIBox2 clipBox,
				   WlzErrorNum *dstErrNum)
{
  int		itvCount,
  		lnCount,
		lnIdx;
  WlzObject	*dstObj = NULL;
  WlzDomain	dstDom,
  		srcDom;
  WlzInterval	*dstItv0,
  		*dstItv1,
		*dstItv2,
  		*srcItv;
  WlzIntervalLine *itvLn;
  WlzIBox2	relBox,
  		lnBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzClipObjToBox2D FE 0x%lx {%d %d %d %d} 0x%lx\n",
	   (unsigned long )srcObj,
	   clipBox.xMin, clipBox.yMin, clipBox.xMax, clipBox.yMax,
	   (unsigned long )dstErrNum));
  dstDom.core = NULL;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    WLZ_DBG((WLZ_DBG_LVL_2),
	    ("WlzClipObjToBox2D 01 %d\n",
	     (int )(srcObj->type)));
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
        break;
      case WLZ_2D_DOMAINOBJ:
	if((srcDom = srcObj->domain).core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if((srcDom.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
		(srcDom.core->type != WLZ_INTERVALDOMAIN_RECT))
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  clipBox.xMin = WLZ_MAX(clipBox.xMin, srcDom.i->kol1); 
	  clipBox.yMin = WLZ_MAX(clipBox.yMin, srcDom.i->line1); 
	  clipBox.xMax = WLZ_MIN(clipBox.xMax, srcDom.i->lastkl); 
	  clipBox.yMax = WLZ_MIN(clipBox.yMax, srcDom.i->lastln); 
	  if((clipBox.xMin <= clipBox.xMax) && (clipBox.yMin <= clipBox.yMax))
	  {
	    WLZ_DBG((WLZ_DBG_LVL_2),
		    ("WlzClipObjToBox2D 02 %d\n",
		     (int )(srcDom.core->type)));
	    switch(srcDom.core->type)
	    {
	      case WLZ_INTERVALDOMAIN_INTVL:
		itvCount = 0;
		lnCount = clipBox.yMax - clipBox.yMin;
		itvLn = srcDom.i->intvlines + clipBox.yMin - srcDom.i->line1;
		while(lnCount-- >= 0)
		{
		  itvCount += itvLn->nintvs;
		  ++itvLn;
		}
		if(((dstItv0 = (WlzInterval *)AlcMalloc(sizeof(WlzInterval) *
					 (unsigned long )itvCount)) != NULL) &&
		   ((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
						      clipBox.yMin,
						      clipBox.yMax,
						      clipBox.xMin,
						      clipBox.xMax,
						      &errNum)) != NULL))
		{
		  dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr,
		  				       (void *)dstItv0, NULL);
		  itvLn = srcDom.i->intvlines + clipBox.yMin - srcDom.i->line1;
		  relBox = clipBox;
		  relBox.xMin -= srcDom.i->kol1;
		  relBox.xMax -= srcDom.i->kol1;
		  lnBox.xMin = srcDom.i->lastkl;
		  lnBox.xMax = 0;
		  dstItv1 = dstItv0;
		  lnIdx = clipBox.yMin;
		  lnCount = clipBox.yMax - clipBox.yMin + 1;
		  while(lnCount-- > 0)
		  {
		    itvCount = itvLn->nintvs;
		    srcItv = itvLn->intvs;
		    dstItv2 = dstItv1;
		    while(itvCount-- > 0)
		    {
		      if(srcItv->ileft <= relBox.xMax)
		      {
			if(srcItv->iright >= relBox.xMin)
			{
			  dstItv1->ileft = WLZ_MAX(relBox.xMin, srcItv->ileft);
			  dstItv1->iright = WLZ_MIN(relBox.xMax,
			  			    srcItv->iright);
			  lnBox.xMin = WLZ_MIN(lnBox.xMin, dstItv1->ileft);
			  lnBox.xMax = WLZ_MAX(lnBox.xMax, dstItv1->iright);
			  dstItv1++;
			}
			srcItv++;
		      }
		    }
		    WlzMakeInterval(lnIdx, dstDom.i, (int )(dstItv1 - dstItv2),
				    dstItv2);
		    ++itvLn;
		    ++lnIdx;
		  } 
		  dstDom.i->lastkl = srcDom.i->kol1 + lnBox.xMax;
		  dstDom.i->kol1 = lnBox.xMin + srcDom.i->kol1;
		  dstItv2 = dstItv0;
		  itvCount = (int )(dstItv1 - dstItv0);
		  while(itvCount-- > 0)
		  {
		    dstItv2->ileft -= lnBox.xMin;
		    dstItv2->iright -= lnBox.xMin;
		    dstItv2++;
		  }
		}
		break;
	      case WLZ_INTERVALDOMAIN_RECT:
		dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					         clipBox.yMin, clipBox.yMax,
					         clipBox.xMin, clipBox.xMax,
						 &errNum);
		break;
	      default:
	        errNum = WLZ_ERR_DOMAIN_TYPE;
		break;
	    }
	    if(dstDom.i)
	    {
	      dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom,
				   srcObj->values, NULL, NULL, &errNum);
	    }
	  }
	  else  /* No intersection between the clip box and the given domain */
	  {
	    dstObj = WlzMakeEmpty(&errNum);
	  }
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      WlzFreeObj(dstObj);
      dstObj = NULL;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_2),
	  ("WlzClipObjToBox2D 03 %d\n",
	   (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzClipObjToBox2D FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/*!
* \return	New object with clipped domain or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Clips the given object's domain so that it does not lie outside
* 		the given 3D axis alaigned clip box.
* \param	srcObj			Given source object.
* \param	clipBox			Clip box.
* \param	dstErrNum		Destination pointer for error, may be
* 					NULL.
*/
WlzObject	*WlzClipObjToBox3D(WlzObject *srcObj, WlzIBox3 clipBox,
				   WlzErrorNum *dstErrNum)
{
  int		dstPlaneIdx,
		srcPlaneIdx,
		planeCount;
		WlzObject	*dstObj = NULL,
		*srcObj2D = NULL,
		*dstObj2D = NULL;
  WlzDomain	dstDom,
		srcDom,
		srcDom2D;
  WlzValues	dstValues,
		srcValues2D;
  WlzIBox2	clipBox2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzClipObjToBox3D FE 0x%lx {%d %d %d %d %d %d} 0x%lx\n",
	   (unsigned long )srcObj,
	   clipBox.xMin, clipBox.yMin, clipBox.xMax, clipBox.yMax,
	   clipBox.zMin, clipBox.zMax,
	   (unsigned long )dstErrNum));
  dstDom.core = NULL;
  dstValues.core = NULL;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    WLZ_DBG((WLZ_DBG_LVL_2),
	    ("WlzClipObjToBox3D 01 %d\n",
	     (int )(srcObj->type)));
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ:
	clipBox2D.xMin = clipBox.xMin;
	clipBox2D.xMax = clipBox.xMax;
	clipBox2D.yMin = clipBox.yMin;
	clipBox2D.yMax = clipBox.yMax;
	dstObj = WlzClipObjToBox2D(srcObj, clipBox2D, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
	if((srcDom = srcObj->domain).core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  clipBox.xMin = WLZ_MAX(clipBox.xMin, srcDom.p->kol1); 
	  clipBox.yMin = WLZ_MAX(clipBox.yMin, srcDom.p->line1); 
	  clipBox.zMin = WLZ_MAX(clipBox.zMin, srcDom.p->plane1); 
	  clipBox.xMax = WLZ_MIN(clipBox.xMax, srcDom.p->lastkl); 
	  clipBox.yMax = WLZ_MIN(clipBox.yMax, srcDom.p->lastln); 
	  clipBox.zMax = WLZ_MIN(clipBox.zMax, srcDom.p->lastpl); 
	  if((clipBox.xMin <= clipBox.xMax) &&
	     (clipBox.yMin <= clipBox.yMax) &&
	     (clipBox.zMin <= clipBox.zMax))
	  {
	    srcDom2D.core = NULL;
	    srcValues2D.core = NULL;
	    clipBox2D.xMin = clipBox.xMin;
	    clipBox2D.xMax = clipBox.xMax;
	    clipBox2D.yMin = clipBox.yMin;
	    clipBox2D.yMax = clipBox.yMax;
	    dstDom.p =  WlzMakePlaneDomain(srcDom.p->type,
					   clipBox.zMin, clipBox.zMax,
					   clipBox.yMin, clipBox.yMax,
					   clipBox.xMin,
					   clipBox.xMax, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(srcObj->values.core == NULL) /* For backwards compatibility */
	      {
		dstValues.vox = NULL;
	      }
	      else
	      {
		if(srcObj->values.core->type == WLZ_EMPTY_OBJ)
		{
		  dstValues.vox = (WlzVoxelValues *)WlzMakeEmpty(&errNum);
		}
		else
		{
		  dstValues.vox = WlzMakeVoxelValueTb(srcObj->values.vox->type,
						      clipBox.zMin,
						      clipBox.zMax,
						      WlzGetBackground(srcObj,
								       NULL),
						      NULL, &errNum);
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, srcDom2D,
				     srcValues2D, NULL, NULL, &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dstPlaneIdx = 0;
	      srcPlaneIdx = clipBox.zMin - srcDom.p->plane1;
	      planeCount = clipBox.zMax - clipBox.zMin + 1;
 	      if((dstValues.core == NULL) ||
	         (dstValues.core->type == WLZ_EMPTY_OBJ))
	      {
		while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
		{
		  srcObj2D->domain = *(srcObj->domain.p->domains +
		  		       srcPlaneIdx);
		  if(srcObj2D->domain.core == NULL)      /* Backwards compat */
		  {
		    *(dstDom.p->domains + dstPlaneIdx) =
					     WlzAssignDomain(srcObj2D->domain,
							     NULL);
		  }
		  else
		  {
		    if(((dstObj2D = WlzClipObjToBox2D(srcObj2D, clipBox2D,
						      &errNum)) != NULL) &&
		       (errNum == WLZ_ERR_NONE))
		    {
		      *(dstDom.p->domains + dstPlaneIdx) =
					     WlzAssignDomain(dstObj2D->domain,
							     NULL);
		      WlzFreeObj(dstObj2D);
		    }
		  }
		  ++srcPlaneIdx;
		  ++dstPlaneIdx;
		}
	      }
	      else
	      {
		while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
		{
		  srcObj2D->domain = *(srcObj->domain.p->domains + srcPlaneIdx);
		  srcObj2D->values = *(srcObj->values.vox->values +
				       srcPlaneIdx);
		  if(srcObj2D->domain.core == NULL)      /* Backwards compat */
		  {
		    *(dstDom.p->domains + dstPlaneIdx) =
					     WlzAssignDomain(srcObj2D->domain,
							     NULL);
		    *(dstValues.vox->values + dstPlaneIdx) =
					     WlzAssignValues(srcObj2D->values,
							     NULL);
		  }
		  else
		  {
		    if(((dstObj2D = WlzClipObjToBox2D(srcObj2D, clipBox2D,
						      &errNum)) != NULL) &&
		       (errNum == WLZ_ERR_NONE))
		    {
		      *(dstDom.p->domains + dstPlaneIdx) =
					     WlzAssignDomain(dstObj2D->domain,
							     NULL);
		      *(dstValues.vox->values + dstPlaneIdx) =
					     WlzAssignValues(dstObj2D->values,
							     NULL);
		      WlzFreeObj(dstObj2D);
		    }
		  }
		  ++srcPlaneIdx;
		  ++dstPlaneIdx;
		}
	      }
	      srcObj2D->domain.core = NULL;
	      srcObj2D->values.core = NULL;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dstDom.p->voxel_size[0] = srcObj->domain.p->voxel_size[0];
	      dstDom.p->voxel_size[1] = srcObj->domain.p->voxel_size[1];
	      dstDom.p->voxel_size[2] = srcObj->domain.p->voxel_size[2];
	      WlzStandardPlaneDomain(dstDom.p, dstValues.vox);
	      dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues, NULL,
		  		   NULL, &errNum);

	    }
	  }
	  else
	  {
	    dstObj = WlzMakeEmpty(&errNum);
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(srcObj2D)
  {
    WlzFreeObj(srcObj2D);
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    dstObj = NULL;
  }
  WLZ_DBG((WLZ_DBG_LVL_2),
      ("WlzClipObjToBox3D 03 %d\n",
       (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
      ("WlzClipObjToBox3D FX 0x%lx\n",
       (unsigned long )dstObj));
  return(dstObj);
}
