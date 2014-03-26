#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzIterate_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzIterate.c
* \author       Bill Hill
* \date         April 2011
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
* \brief	Functions for iteration through Woolz objects.
* \ingroup	WlzDomainOps
*/

#include <string.h>
#include <Wlz.h>

static WlzIterateWSpace 	*WlzIterateMakeWSp();
static WlzErrorNum 		WlzIterateInitDomObj2D(
				  WlzIterateWSpace *itWSp,
				  WlzObject *obj,
				  WlzRasterDir dir,
				  int grey);
static WlzErrorNum 		WlzIterateInitDomObj3D(
				  WlzIterateWSpace *itWSp,
				  WlzObject *obj,
				  WlzRasterDir dir,
				  int grey);
static WlzErrorNum 		WlzIterateDomObj2D(
				  WlzIterateWSpace *itWSp);
static WlzErrorNum 		WlzIterateDomObj3D(
				  WlzIterateWSpace *itWSp);

/*!
* \return	Returns a 2D raster direction.
* \ingroup	WlzDomainOps
* \brief	Masks a 2 or 3D raster direction so that it is valid for 2D.
* \param	dir			Given raster direction.
*/
WlzRasterDir	WlzRasterDir2D(WlzRasterDir dir)
{
  dir = dir & WLZ_RASTERDIR_DLDC;
  return(dir);
}

/*!
* \return	A new domain object iteration workspace.
* \ingroup	WlzDomainOps
* \brief	Creates a new domain object iteration workspace that may be
* 		used to iterate through either a 2 or 3D domain object.
* \param	obj			Given object which must be either of
* 					the type WLZ_2D_DOMAINOBJ or
* 					WLZ_3D_DOMAINOBJ.
* \param	dir			Raster scan direction.
* \param	grey			Non-zero if grey value access is needed
* 					(not required if the object does not
* 					have grey values).
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIterateWSpace *WlzIterateInit(WlzObject *obj, WlzRasterDir dir,
                                 int grey, WlzErrorNum *dstErr)
{
  WlzIterateWSpace *itWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((itWSp = WlzIterateMakeWSp()) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	dir = WlzRasterDir2D(dir);
	errNum = WlzIterateInitDomObj2D(itWSp, obj, dir, grey);
	break;
      case WLZ_3D_DOMAINOBJ:
	errNum = WlzIterateInitDomObj3D(itWSp, obj, dir, grey);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzIterateWSpFree(itWSp);
    itWSp = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(itWSp);
}

/*!
* \ingroup	WlzDomainOps
* \brief	Frees the given iteration workspace including it's
* 		interval and grey workspaces. The object and 2D object
* 		are freed if non NULL.
* \param	itWSp			Given iteration workspace.
*/
void		WlzIterateWSpFree(WlzIterateWSpace *itWSp)
{
  if(itWSp)
  {
    if(itWSp->iWSp && itWSp->gWSp)
    {
      WlzEndGreyScan(itWSp, itWSp->gWSp);
    }
    AlcFree(itWSp->iWSp);
    AlcFree(itWSp->gWSp);
    (void )WlzFreeObj(itWSp->obj2D);
    AlcFree(itWSp);
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Iterates through an object for which the given workspace was
* 		initialised. This function will return WLZ_ERR_EOO when there
* 		are no more pixels/voxels remaining.
* \param	itWSp				Given iteration workspace.
*/
WlzErrorNum	WlzIterate(WlzIterateWSpace *itWSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(itWSp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(itWSp->obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(itWSp->obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	errNum = WlzIterateDomObj2D(itWSp);
	break;
      case WLZ_3D_DOMAINOBJ:
	errNum = WlzIterateDomObj3D(itWSp);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	New iteration workspace data structure.
* \ingroup	WlzDomainOps
* \brief	Just alocates a new iteration workspace data structure.
*/
static WlzIterateWSpace *WlzIterateMakeWSp()
{
  WlzIterateWSpace *itWSp = NULL;

  itWSp = AlcMalloc(sizeof(WlzIterateWSpace));
  return(itWSp);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Initialises an iteration workspace data structure for 2D
* 		domain objects.
* \param	itWSp			Iteration workspace data structure
* 					to be initialised.
* \param	obj			Given object for which to initialise
* 					the iteration workspace.
* \param	dir			Raster scan direction.
* \param	grey			Non-zero if grey value access is needed
* 					(not required if the object does not
* 					have grey values).
*/
static WlzErrorNum WlzIterateInitDomObj2D(WlzIterateWSpace *itWSp,
				WlzObject *obj, WlzRasterDir dir, int grey)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  (void )memset(itWSp, 0, sizeof(WlzIterateWSpace));
  itWSp->dir = dir;
  itWSp->obj = obj;
  if((itWSp->iWSp = AlcMalloc(sizeof(WlzIntervalWSpace))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )memset(itWSp->iWSp, 0, sizeof(WlzIntervalWSpace));
    itWSp->gType = WLZ_GREY_ERROR;
    if(grey && obj->values.core)
    {
      itWSp->grey = grey;
      if((itWSp->gWSp = AlcMalloc(sizeof(WlzGreyWSpace))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        (void )memset(itWSp->gWSp, 0, sizeof(WlzGreyWSpace));
	errNum = WlzInitGreyRasterScan(itWSp->obj, itWSp->iWSp, itWSp->gWSp,
	                               dir, 0);
      }
    }
    else
    {
      errNum = WlzInitRasterScan(itWSp->obj, itWSp->iWSp, dir);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Initialises an iteration workspace data structure for 2D
* 		domain objects.
* \param	itWSp			Iteration workspace data structure
* 					to be initialised.
* \param	obj			Given object for which to initialise
* 					the iteration workspace.
* \param	dir			Raster scan direction.
* \param	grey			Non-zero if grey value access is needed
* 					(not required if the object does not
* 					have grey values).
*/
static WlzErrorNum WlzIterateInitDomObj3D(WlzIterateWSpace *itWSp,
				WlzObject *obj, WlzRasterDir dir, int grey)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  (void )memset(itWSp, 0, sizeof(WlzIterateWSpace));
  itWSp->dir = dir;
  itWSp->obj = obj;
  if((itWSp->iWSp = AlcCalloc(1, sizeof(WlzIntervalWSpace))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    itWSp->gType = WLZ_GREY_ERROR;
    itWSp->iWSp->linrmn = -1;
    if(grey && obj->values.core)
    {
      itWSp->grey = grey;
      if((itWSp->gWSp = AlcCalloc(1, sizeof(WlzGreyWSpace))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Iterates through a 2D domain object.
* \param	itWSp			Given iteration workspace.
*/
static WlzErrorNum WlzIterateDomObj2D(WlzIterateWSpace *itWSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(itWSp->itvPos);
  if((itWSp->iWSp->colrmn - itWSp->itvPos) <= 0)
  {
    if(itWSp->grey)
    {
      errNum = WlzNextGreyInterval(itWSp->iWSp);
    }
    else
    {
      errNum = WlzNextInterval(itWSp->iWSp);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      itWSp->itvPos = 0;
      itWSp->pos.vtX = itWSp->iWSp->lftpos;
      itWSp->pos.vtY = itWSp->iWSp->linpos;
      if(itWSp->grey)
      {
        itWSp->gType = itWSp->gWSp->pixeltype;
	itWSp->gP = itWSp->gWSp->u_grintptr;
      }
    }
  }
  else
  {
    ++(itWSp->pos.vtX);
    switch(itWSp->gType)
    {
      case WLZ_GREY_LONG:
        ++(itWSp->gP.lnp);
	break;
      case WLZ_GREY_INT:
        ++(itWSp->gP.inp);
	break;
      case WLZ_GREY_SHORT:
        ++(itWSp->gP.shp);
	break;
      case WLZ_GREY_UBYTE:
        ++(itWSp->gP.ubp);
	break;
      case WLZ_GREY_FLOAT:
        ++(itWSp->gP.flp);
	break;
      case WLZ_GREY_DOUBLE:
        ++(itWSp->gP.dbp);
	break;
      case WLZ_GREY_BIT:
        ++(itWSp->gP.ubp);
	break;
      case WLZ_GREY_RGBA:
        ++(itWSp->gP.rgbp);
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Iterates through a 3D domain object.
* \param	itWSp			Given iteration workspace.
*/
static WlzErrorNum WlzIterateDomObj3D(WlzIterateWSpace *itWSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(itWSp->itvPos);
  if((itWSp->iWSp == NULL) ||
     (itWSp->iWSp->colrmn - itWSp->itvPos) <= 0)
  {
    if((itWSp->iWSp == NULL) ||
       (itWSp->iWSp->objaddr == NULL))
    {
      errNum = WLZ_ERR_EOO;
    }
    else
    {
      errNum = (itWSp->grey)? WlzNextGreyInterval(itWSp->iWSp):
			      WlzNextInterval(itWSp->iWSp);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      itWSp->itvPos = 0;
      itWSp->pos.vtX = itWSp->iWSp->lftpos;
      itWSp->pos.vtY = itWSp->iWSp->linpos;
      if(itWSp->grey)
      {
	itWSp->gType = itWSp->gWSp->pixeltype;
	itWSp->gP = itWSp->gWSp->u_grintptr;
      }
    }
    else if(errNum == WLZ_ERR_EOO)
    {
      int 	  pCnt,
		  pIdx,
		  pInc;
      WlzPlaneDomain *pDom;

      errNum = WLZ_ERR_NONE;
      pDom = itWSp->obj->domain.p;
      pCnt = pDom->lastpl - pDom->plane1 + 1;
      pInc = (((itWSp->dir) & WLZ_RASTERDIR_DP) == 0)? 1: -1;
      --(itWSp->plnRmn);
      if(itWSp->plnRmn == 0)
      {
	/* No more planes to iterate through. */
	errNum = WLZ_ERR_EOO;
      }
      else if(itWSp->plnRmn < 0)
      {
	/* This is the first plane coming up. */
	itWSp->plnIdx = (pInc > 0)? 0: pCnt - 1;
	itWSp->plnRmn = pCnt;
      }
      else
      {
	itWSp->plnIdx += pInc;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	pIdx = itWSp->plnIdx;
	while(((pInc > 0) && (pIdx < pCnt)) ||
	      ((pInc < 0) && (pIdx >= 0)))
	{
	  WlzDomain	dom;

	  itWSp->plnIdx = pIdx;
	  dom = *(pDom->domains + pIdx);
	  if((dom.core != NULL) && (dom.core->type != WLZ_EMPTY_DOMAIN))
	  {
	    WlzValues val;

	    val.core = NULL;
	    itWSp->pos.vtZ = pDom->plane1 + pIdx;
	    if(itWSp->grey)
	    {
	      val = *(itWSp->obj->values.vox->values + pIdx);
	    }
	    (void )WlzFreeObj(itWSp->obj2D);
	    itWSp->obj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom, val,
				       NULL, NULL, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      WlzRasterDir dir2D;

              itWSp->itvPos = 0;
	      dir2D = WlzRasterDir2D(itWSp->dir);
	      if(itWSp->grey)
	      {
                (void )memset(itWSp->gWSp, 0, sizeof(WlzGreyWSpace));
		errNum = WlzInitGreyRasterScan(itWSp->obj2D,
					       itWSp->iWSp, itWSp->gWSp,
					       dir2D, 0);
		{
                  errNum = WlzNextGreyInterval(itWSp->iWSp);
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  itWSp->gType = itWSp->gWSp->pixeltype;
		  itWSp->gP = itWSp->gWSp->u_grintptr;
		}
	      }
	      else
	      {
		errNum = WlzInitRasterScan(itWSp->obj2D, itWSp->iWSp, dir2D);
		if(errNum == WLZ_ERR_NONE)
		{
                  errNum = WlzNextInterval(itWSp->iWSp);
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      itWSp->pos.vtX = itWSp->iWSp->lftpos;
	      itWSp->pos.vtY = itWSp->iWSp->linpos;
	    }
	    break;
	  }
	  pIdx += pInc;
	}
	itWSp->itvPos = 0;
	itWSp->plnRmn = (pInc > 0)? pCnt - itWSp->plnIdx: itWSp->plnIdx + 1;
      }
    }
  }
  else
  {
    ++(itWSp->pos.vtX);
    if(itWSp->grey)
    {
      switch(itWSp->gType)
      {
	case WLZ_GREY_LONG:
	  ++(itWSp->gP.lnp);
	  break;
	case WLZ_GREY_INT:
	  ++(itWSp->gP.inp);
	  break;
	case WLZ_GREY_SHORT:
	  ++(itWSp->gP.shp);
	  break;
	case WLZ_GREY_UBYTE:
	  ++(itWSp->gP.ubp);
	  break;
	case WLZ_GREY_FLOAT:
	  ++(itWSp->gP.flp);
	  break;
	case WLZ_GREY_DOUBLE:
	  ++(itWSp->gP.dbp);
	  break;
	case WLZ_GREY_BIT:
	  ++(itWSp->gP.ubp);
	  break;
	case WLZ_GREY_RGBA:
	  ++(itWSp->gP.rgbp);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}
