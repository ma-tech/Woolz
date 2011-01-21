#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreyValue_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzGreyValue.c
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
* \brief	Provides functions for random access to the grey values
* 		of 2D and 3D domain objects.
* \ingroup	WlzAccess
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <limits.h>
#include <Wlz.h>

static void			WlzGreyValueSetBkdP(
				  WlzGreyV *gVP,
				  WlzGreyP *gPP,
				  WlzGreyType gType,
				  WlzGreyV val);
static void			WlzGreyValueSetBkdPN(
				  WlzGreyV *gVP,
				  WlzGreyP *gPP,
				  WlzGreyType gType,
				  WlzGreyV val,
				  int count);
static void			WlzGreyValueSetGreyP(
				  WlzGreyV *gVP,
				  WlzGreyP *gPP,
				  WlzGreyType gType,
	  			  WlzGreyP baseGVP,
				  size_t offset);
static void			WlzGreyValueComputeGreyP2D(
				  WlzGreyP *baseGVP,
				  size_t *offset,
				  WlzGreyValueWSpace *gVWSp,
				  int line,
				  int kol);
static void			WlzGreyValueComputeGreyPTiled2D(
				  WlzGreyP *baseGVP,
				  size_t *offset,
				  WlzTiledValues *tVal,
				  int line,
				  int kol);
static void			WlzGreyValueComputeGreyPTiled3D(
				  WlzGreyP *baseGVP,
				  size_t *offset,
				  WlzTiledValues *tVal,
				  int plane,
				  int line,
				  int kol);
static void			WlzGreyValueGet2D1(
				  WlzGreyValueWSpace *gVWSp,
				  int line,
				  int kol);
static void			WlzGreyValueGet3D1(
				  WlzGreyValueWSpace *gVWSp,
				  int plane,
				  int line,
				  int kol);
static void			WlzGreyValueGet3DTiled(
				  WlzGreyValueWSpace *gVWSp,
				  int plane,
				  int line,
				  int kol);
static void			WlzGreyValueGet2DCon(
				  WlzGreyValueWSpace *gVWSp,
				  int line,
				  int kol);
static void			WlzGreyValueGet3DCon(
				  WlzGreyValueWSpace *gVWSp,
				  int plane,
				  int line,
				  int kol);
static void			WlzGreyValueGet3DConTiled(
				  WlzGreyValueWSpace *gVWSp,
				  int plane,
				  int line,
				  int kol);
static void			WlzGreyValueGetTransCon(
				  WlzGreyValueWSpace *gVWSp,
				  int plane,
				  int line,
				  int kol);
/*!
* \return	Grey value work space or NULL on error.
* \ingroup	WlzAccess
* \brief	Creates a grey value work space from the given object.
*		The resulting grey value work space should be freed
*		using WlzGreyValueFreeWSp().
* \param	obj			Given object.
* \param	dstErrNum		Destination error pointer, may be NULL.
*/
WlzGreyValueWSpace *WlzGreyValueMakeWSp(WlzObject *obj,
					WlzErrorNum *dstErrNum)
{
  int		planeIdx, numPlanes;
  WlzDomain	*planeDomains;
  WlzValues	*planeValues;
  WlzObjectType	gTabType0,
	        gTabType1 = WLZ_DUMMY_ENTRY;
  WlzGreyType	gType0,
		gType1 = WLZ_GREY_ERROR;
  WlzPixelV	bkdPix;
  WlzAffineTransform *trans0,
  		*trans1;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
  	   ("WlzGreyValueMakeWSp FE %p %p\n",
	    obj, dstErrNum));
  if(obj == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(obj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else if((gVWSp = (WlzGreyValueWSpace *)AlcCalloc(1,
					  sizeof(WlzGreyValueWSpace))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  gVWSp->objType = obj->type;
	  gVWSp->values = obj->values;
	  gVWSp->gType = WlzGreyTableTypeToGreyType(obj->values.core->type,
	  					    &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    gVWSp->gTabType = WlzGreyTableTypeToTableType(
	    				obj->values.core->type,
					&errNum);
	    gVWSp->gTabType2D = gVWSp->gTabType;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    switch(gVWSp->gTabType2D)
	    {
	      case WLZ_GREY_TAB_RAGR:
	        gVWSp->values2D = obj->values;
		bkdPix = gVWSp->values2D.v->bckgrnd;
		break;
	      case WLZ_GREY_TAB_RECT:
	        gVWSp->values2D = obj->values;
		bkdPix = gVWSp->values2D.r->bckgrnd;
		break;
	      case WLZ_GREY_TAB_INTL:
	        gVWSp->values2D = obj->values;
		bkdPix = gVWSp->values2D.i->bckgrnd;
		break;
	      case WLZ_GREY_TAB_TILED:
	        gVWSp->values2D = obj->values;
		bkdPix = gVWSp->values2D.t->bckgrnd;
		break;
	      default:
		errNum = WLZ_ERR_VALUES_TYPE;
		break;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    switch(gVWSp->gType)
	    {
	      case WLZ_GREY_LONG:
	      case WLZ_GREY_INT:
	      case WLZ_GREY_SHORT:
	      case WLZ_GREY_UBYTE:
	      case WLZ_GREY_FLOAT:
	      case WLZ_GREY_DOUBLE:
	      case WLZ_GREY_RGBA:
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	  }
	  if((gVWSp->gType != bkdPix.type) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WLZ_ERR_VALUES_DATA;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  gVWSp->gBkd = bkdPix.v;
	  gVWSp->domain = obj->domain;
	  gVWSp->iDom2D = obj->domain.i;
	  if((gVWSp->iDom2D->type != WLZ_INTERVALDOMAIN_INTVL) &&
	     (gVWSp->iDom2D->type != WLZ_INTERVALDOMAIN_RECT))
	  {
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	  }
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(obj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else if(obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if((gVWSp = (WlzGreyValueWSpace *)AlcCalloc(1,
					  sizeof(WlzGreyValueWSpace))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  gVWSp->objType = obj->type;
	  gVWSp->domain = obj->domain;
	  gVWSp->values = obj->values;
	  if(WlzGreyTableIsTiled(obj->values.core->type) == WLZ_GREY_TAB_TILED)
	  {
	    gVWSp->gTabType = WLZ_GREY_TAB_TILED;
	    gVWSp->gTabType2D = WLZ_GREY_TAB_TILED;
	    gVWSp->gType = WlzGreyTableTypeToGreyType(obj->values.core->type,
		                                      NULL);
	    gVWSp->plane = obj->domain.p->plane1;
	    gVWSp->iDom2D = (*(obj->domain.p->domains)).i;
	    gVWSp->gBkd = obj->values.t->bckgrnd.v;
	    if(obj->values.t->bckgrnd.type != gVWSp->gType)
	    {
	      errNum = WLZ_ERR_VALUES_DATA;
	    }
	  }
	  else
	  {
	    gVWSp->gTabType = WLZ_VOXELVALUETABLE_GREY;
	    /* Put in a list of grey-table types - purely for efficiency
	       to avoid re-computing the type for each voxel request.
	       When table types are no longer computed this in principle
	       could go */
	    numPlanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
	    if((gVWSp->gTabTypes3D = (WlzObjectType *)
			  AlcCalloc(numPlanes, sizeof(WlzObjectType))) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      WlzObjectType	*gTabTypes3D;

	      gVWSp->values2D.core = NULL;
	      planeIdx = obj->domain.p->plane1;
	      planeDomains = obj->domain.p->domains;
	      planeValues = obj->values.vox->values;
	      gTabTypes3D = gVWSp->gTabTypes3D;
	      while((planeIdx <= obj->domain.p->lastpl) &&
		    (errNum == WLZ_ERR_NONE))
	      {
		if(((*planeValues).core) &&
		   ((*planeValues).core->type != WLZ_EMPTY_OBJ))
		{
		  gType0 = WlzGreyTableTypeToGreyType((*planeValues).core->type,
		                                      &errNum);
		  if(errNum == WLZ_ERR_NONE)
		  {
		    gTabType0 = WlzGreyTableTypeToTableType((*planeValues).
			                                    core->type,
							    &errNum);
		    *gTabTypes3D = gTabType0;
		  }
		  if(errNum == WLZ_ERR_NONE)
		  {
		    if((gType0 != gType1) && (gType1 != WLZ_GREY_ERROR))
		    {
		      errNum = WLZ_ERR_VALUES_DATA;
		    }
		    else
		    {
		      gType1 = gType0;
		      gTabType1 = gTabType0;
		      if(gVWSp->values2D.core == NULL)
		      {
			gVWSp->plane = planeIdx;
			gVWSp->iDom2D = (*planeDomains).i;
			gVWSp->values2D = *planeValues;
		      }
		    }
		  }
		}
		++planeIdx;
		++planeDomains;
		++planeValues;
		++gTabTypes3D;
	      }
	      if((gType1 == WLZ_GREY_ERROR) || (gTabType1 == WLZ_DUMMY_ENTRY))
	      {
		errNum = WLZ_ERR_VALUES_DATA;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      gVWSp->gType = gType0;
	      gVWSp->gTabType2D = gVWSp->gTabTypes3D[0];
	      switch(gVWSp->gTabType2D)
	      {
		case WLZ_GREY_TAB_RAGR:
		  bkdPix = gVWSp->values2D.v->bckgrnd;
		  break;
		case WLZ_GREY_TAB_RECT:
		  bkdPix = gVWSp->values2D.r->bckgrnd;
		  break;
		case WLZ_GREY_TAB_INTL:
		  bkdPix = gVWSp->values2D.i->bckgrnd;
		  break;
		case WLZ_GREY_TAB_TILED:
		  bkdPix = gVWSp->values2D.t->bckgrnd;
		  break;
		default:
		  errNum = WLZ_ERR_VALUES_TYPE;
		  break;
	      }
	      if((gVWSp->gType != bkdPix.type) && (errNum == WLZ_ERR_NONE))
	      {
		errNum = WLZ_ERR_VALUES_DATA;
	      }
	      else if((gVWSp->iDom2D->type != WLZ_INTERVALDOMAIN_INTVL) &&
		      (gVWSp->iDom2D->type != WLZ_INTERVALDOMAIN_RECT))
	      {
		errNum = WLZ_ERR_DOMAIN_TYPE;
	      }
	      else
	      {
		gVWSp->gBkd = bkdPix.v;
	      }
	    }
	  }
	}
	break;
      case WLZ_TRANS_OBJ:
	trans0 = NULL;
	while((errNum == WLZ_ERR_NONE) && (obj->type == WLZ_TRANS_OBJ))
	{
	  if(trans0 == NULL)
	  {
	    trans0 = WlzAffineTransformCopy(obj->domain.t, &errNum);
	  }
	  else
	  {
	    if((trans1 = WlzAffineTransformProduct(trans0, obj->domain.t,
	    				           &errNum)) != NULL)
	    {
	      WlzFreeAffineTransform(trans0);
	      trans0 = trans1;
	      trans1 = NULL;
	    }
	  }
	  if((obj = obj->values.obj) == NULL)
	  {
	    errNum = WLZ_ERR_OBJECT_NULL;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  gVWSp->invTrans = WlzAffineTransformInverse(trans0, &errNum);
	}
	if(trans0)
	{
	  WlzFreeAffineTransform(trans0);
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if((errNum != WLZ_ERR_NONE) && (gVWSp != NULL))
  {
    WlzGreyValueFreeWSp(gVWSp);
    gVWSp = NULL;
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzGreyValueMakeWSp FX %p\n",
	   gVWSp));
  return(gVWSp);
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Free's the given grey value work space created by
* 		WlzGreyValueMakeWSp().
* \param	gVWSp			Given grey value work space.
*/
void		WlzGreyValueFreeWSp(WlzGreyValueWSpace *gVWSp)
{
  WLZ_DBG((WLZ_DBG_LVL_1),
  	   ("WlzGreyValueFreeWSp FE %p\n",
	    gVWSp));
  if(gVWSp)
  {
    (void )WlzFreeAffineTransform(gVWSp->invTrans);
    AlcFree((void *)(gVWSp->gTabTypes3D));
    AlcFree(gVWSp);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzGreyValueFreeWSp FX\n"));
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets a single grey value/pointer for the given point
*		from the object with which the given work space was
*               initialised.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane (z) coordinate of point.
* \param	line			Line (y) coordinate of point.
* \param	kol			Column (x) coordinate of point.
*/
void		WlzGreyValueGet(WlzGreyValueWSpace *gVWSp,
			        double plane, double line, double kol)
{
  if(gVWSp)
  {
    switch(gVWSp->objType)
    {
      case WLZ_2D_DOMAINOBJ:
	if(gVWSp->invTrans)
	{
	  WlzDVertex2	vtx2;

	  vtx2.vtX = kol;
	  vtx2.vtY = line;
	  vtx2 = WlzAffineTransformVertexD2(gVWSp->invTrans, vtx2, NULL);
	  kol = vtx2.vtX;
	  line = vtx2.vtY;
	}
	WlzGreyValueGet2D1(gVWSp, WLZ_NINT(line), WLZ_NINT(kol));
	break;
      case WLZ_3D_DOMAINOBJ:
	if(gVWSp->invTrans)
	{
	  WlzDVertex3	vtx3;

	  vtx3.vtX = kol;
	  vtx3.vtY = line;
	  vtx3.vtZ = plane;
	  vtx3 = WlzAffineTransformVertexD3(gVWSp->invTrans, vtx3, NULL);
	  kol = vtx3.vtX;
	  line = vtx3.vtY;
	  plane = vtx3.vtZ;
	}
	if(gVWSp->gTabType == WLZ_GREY_TAB_TILED)
	{
	  WlzGreyValueGet3DTiled(gVWSp,
			     WLZ_NINT(plane), WLZ_NINT(line), WLZ_NINT(kol));
	}
	else
	{
	  WlzGreyValueGet3D1(gVWSp,
			     WLZ_NINT(plane), WLZ_NINT(line), WLZ_NINT(kol));
	}
	break;
      default:
	break;
    }
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets the four/eight connected grey values/pointers for
*               the given point which lie at:
*		\verbatim
                  (p, k, l),          (p, k + 1, l),
                  (p, k, l + 1),      (p, k + 1, l + 1),
                  (p + 1, k, l),      (p + 1, k + 1, l),
                  (p + 1, k, l + 1),  (p + 1, k + 1, l + 1).
*		\endverbatim
* \param	gVWSp			Grey value work space.
* \param	plane			Plane (z) coordinate of point.
* \param	line			Line (y) coordinate of point.
* \param	kol			Column (x) coordinate of point.
*/
void		WlzGreyValueGetCon(WlzGreyValueWSpace *gVWSp,
			           double plane, double line, double kol)
{
  if(gVWSp)
  {
    if(gVWSp->invTrans)
    {
      WlzGreyValueGetTransCon(gVWSp, (int )plane, (int )line, (int )kol);
    }
    else
    {
      /* I don't understand what this can be for?
      * I've used kol, line and plane instead of these below. BH 20100326. 
      *
      ikol = (int) ((kol < 0.0) ? kol - 1.0 : kol);
      iline = (int) ((kol < 0.0) ? line - 1.0 : line);
      iplane = (int) ((kol < 0.0) ? plane - 1.0 : plane);
      */
      switch(gVWSp->objType)
      {
	case WLZ_2D_DOMAINOBJ:
	  WlzGreyValueGet2DCon(gVWSp, line, kol);
	  break;
	case WLZ_3D_DOMAINOBJ:
	  if(gVWSp->gTabType == WLZ_GREY_TAB_TILED)
	  {
	    WlzGreyValueGet3DConTiled(gVWSp, plane, line, kol);
	  }
	  else
	  {
	    WlzGreyValueGet3DCon(gVWSp, plane, line, kol);
	  }
	  break;
	default:
	  break;
      }
    }
  }
}

/*!
* \return	The workspace's object's grey type.
* \ingroup	WlzAccess
* \brief	Access function to get the object's grey type.
* \param	gVWSp			Grey value work space.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGreyType	WlzGreyValueGetGreyType(WlzGreyValueWSpace *gVWSp,
					WlzErrorNum *dstErr)
{
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gVWSp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    gType = gVWSp->gType;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(gType);
}

/*!
* \return	The grey value as an int.
* \ingroup	WlzAccess
* \brief	Gets a single grey value for the given point from the
*               object with which the given work space was initialised.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane (z) coordinate of point.
* \param	line			Line (y) coordinate of point.
* \param	kol			Column (x) coordinate of point.
*/
int		WlzGreyValueGetI(WlzGreyValueWSpace *gVWSp,
			         double plane, double line, double kol)
{
  int		val = 0;

  if(gVWSp)
  {
    WlzGreyValueGet(gVWSp, plane, line, kol);
    switch(gVWSp->gType)
    {
      case WLZ_GREY_LONG:
        val = WLZ_CLAMP(gVWSp->gVal[0].lnv, LLONG_MIN, LLONG_MAX);
        break;
      case WLZ_GREY_INT:
        val = gVWSp->gVal[0].inv;
        break;
      case WLZ_GREY_SHORT:
        val = gVWSp->gVal[0].shv;
        break;
      case WLZ_GREY_UBYTE:
        val = gVWSp->gVal[0].ubv;
        break;
      case WLZ_GREY_FLOAT:
        val = WLZ_NINT(gVWSp->gVal[0].flv);
	break;
      case WLZ_GREY_DOUBLE:
        val = WLZ_NINT(gVWSp->gVal[0].dbv);
	break;
      case WLZ_GREY_RGBA:
        val = gVWSp->gVal[0].rgbv;
	break;
      default:
        break;
    }
  }
  return(val);
}

/*!
* \return	The grey value as a double.
* \ingroup	WlzAccess
* \brief	Gets a single grey value for the given point from the
*               object with which the given work space was initialised.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane (z) coordinate of point.
* \param	line			Line (y) coordinate of point.
* \param	kol			Column (x) coordinate of point.
*/
double		WlzGreyValueGetD(WlzGreyValueWSpace *gVWSp,
				 double plane, double line, double kol)
{
  double		val = 0;

  if(gVWSp)
  {
    WlzGreyValueGet(gVWSp, plane, line, kol);
    switch(gVWSp->gType)
    {
      case WLZ_GREY_LONG:
        val = gVWSp->gVal[0].lnv;
        break;
      case WLZ_GREY_INT:
        val = gVWSp->gVal[0].inv;
        break;
      case WLZ_GREY_SHORT:
        val = gVWSp->gVal[0].shv;
        break;
      case WLZ_GREY_UBYTE:
        val = gVWSp->gVal[0].ubv;
        break;
      case WLZ_GREY_FLOAT:
        val = gVWSp->gVal[0].flv;
	break;
      case WLZ_GREY_DOUBLE:
        val = gVWSp->gVal[0].dbv;
	break;
      case WLZ_GREY_RGBA:
        val = gVWSp->gVal[0].rgbv;
	break;
      default:
        break;
    }
  }
  return(val);
}

/*! 
* \ingroup      WlzValuesUtils
* \brief        Gets a single grey value/pointer for the given point
*		from the 3D values and domain in the work space.
*
* \param    gVWSp	grey-value work space
* \param    plane	plane coordinate
* \param    line	line coordinate
* \param    kol	column coordinate
*/
void		WlzGreyValueGetDir(WlzGreyValueWSpace *gVWSp,
				   int plane, int line, int kol)
{
  int		tI0,
  		valSet = 0;
  WlzDomain	*domP;
  WlzValues	*valP;

  if((plane >= gVWSp->domain.p->plane1) &&
     (plane <= gVWSp->domain.p->lastpl))
  {
    if(gVWSp->gTabType == WLZ_GREY_TAB_TILED)
    {
      WlzGreyValueGet(gVWSp, plane, line, kol);
      valSet = 1;
    }
    else
    {
      if(gVWSp->plane == plane)
      {
	WlzGreyValueGet2D1(gVWSp, line, kol);
	valSet = 1;
      }
      else
      {
	tI0 = plane - gVWSp->domain.p->plane1;
	domP = gVWSp->domain.p->domains + tI0;
	valP = gVWSp->values.vox->values + tI0;
	/* check for non-NULL domain and valuetable
	   pointers until empty obj consistently implemented */
	/*if(domP && valP)*/
	if(domP && valP && (*domP).core && (*valP).core)
	{
	  gVWSp->plane = plane;
	  gVWSp->iDom2D = (*domP).i;
	  gVWSp->values2D = (*valP);
	  gVWSp->gTabType2D = gVWSp->gTabTypes3D[tI0];
	  WlzGreyValueGet2D1(gVWSp, line, kol);
	  valSet = 1;
	}
      }
    }
  }
  if(valSet == 0)
  {
    WlzGreyValueSetBkdP(gVWSp->gVal, gVWSp->gPtr, gVWSp->gType, gVWSp->gBkd);
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Sets the WlzGreyValueWSpace grey value and pointer
*               to the background value.
* \param	gVP			Ptr to grey values.
* \param	gPP			Ptr to grey pointers.
* \param	gType			Grey type.
* \param	val			Background grey value to set.
*/
static void	WlzGreyValueSetBkdP(WlzGreyV *gVP, WlzGreyP *gPP,
				    WlzGreyType gType,
				    WlzGreyV val)
{
  switch(gType)
  {
    case WLZ_GREY_LONG:
      (*gVP).lnv = val.lnv;
      (*(gPP)).lnp = &((*gVP).lnv);
      break;
    case WLZ_GREY_INT:
      (*gVP).inv = val.inv;
      (*(gPP)).inp = &((*gVP).inv);
      break;
    case WLZ_GREY_SHORT:
      (*gVP).shv = val.shv;
      (*(gPP)).shp = &((*gVP).shv);
      break;
    case WLZ_GREY_UBYTE:
      (*gVP).ubv = val.ubv;
      (*(gPP)).ubp = &((*gVP).ubv);
      break;
    case WLZ_GREY_FLOAT:
      (*gVP).flv = val.flv;
      (*(gPP)).flp = &((*gVP).flv);
      break;
    case WLZ_GREY_DOUBLE:
      (*gVP).dbv = val.dbv;
      (*(gPP)).dbp = &((*gVP).dbv);
      break;
    case WLZ_GREY_RGBA:
      (*gVP).rgbv = val.rgbv;
      (*(gPP)).rgbp = &((*gVP).rgbv);
      break;
    default:
      break;
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Sets the WlzGreyValueWSpace grey values and pointers
*               to the background value.
* \param	gVP			Ptr to grey values.
* \param	gPP			Ptr to grey pointers.
* \param	gType			Grey type.
* \param	val			Background grey value to set.
* \param	count			Number of values/pointers set.
*/
static void	WlzGreyValueSetBkdPN(WlzGreyV *gVP, WlzGreyP *gPP,
				     WlzGreyType gType,
				     WlzGreyV val, int count)
{
  switch(gType)
  {
    case WLZ_GREY_LONG:
      while(count-- > 0)
      {
        *gVP = val;
	(*(gPP++)).lnp = &((*gVP++).lnv);
      }
      break;
    case WLZ_GREY_INT:
      while(count-- > 0)
      {
        *gVP = val;
	(*(gPP++)).inp = &((*gVP++).inv);
      }
      break;
    case WLZ_GREY_SHORT:
      while(count-- > 0)
      {
        *gVP = val;
	(*(gPP++)).shp = &((*gVP++).shv);
      }
      break;
    case WLZ_GREY_UBYTE:
      while(count-- > 0)
      {
        *gVP = val;
	(*(gPP++)).ubp = &((*gVP++).ubv);
      }
      break;
    case WLZ_GREY_FLOAT:
      while(count-- > 0)
      {
        *gVP = val;
	(*(gPP++)).flp = &((*gVP++).flv);
      }
      break;
    case WLZ_GREY_DOUBLE:
      while(count-- > 0)
      {
        *gVP = val;
	(*(gPP++)).dbp = &((*gVP++).dbv);
      }
      break;
    case WLZ_GREY_RGBA:
      while(count-- > 0)
      {
        *gVP = val;
	(*(gPP++)).rgbp = &((*gVP++).rgbv);
      }
      break;
    default:
      break;
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Sets the WlzGreyValueWSpace grey values and pointers
*               to the value obtained using the given base pointer and
*               offset.
* \param	gVP			Ptr to grey values.
* \param	gPP			Ptr to grey pointers.
* \param	gType			Grey type.
* \param	baseGVP			Base pointer.
* \param	offset			Offset from base pointer.
*/
static void	WlzGreyValueSetGreyP(WlzGreyV *gVP, WlzGreyP *gPP,
				     WlzGreyType gType,
	  			     WlzGreyP baseGVP, size_t offset)
{
  WlzGreyP	gP;

  if(baseGVP.v != NULL)
  {
    switch(gType)
    {
      case WLZ_GREY_LONG:
	gP.lnp = baseGVP.lnp + offset;
	*(gPP) = gP;
	(*(gVP)).lnv = *(gP.lnp);
	break;
      case WLZ_GREY_INT:
	gP.inp = baseGVP.inp + offset;
	*(gPP) = gP;
	(*(gVP)).inv = *(gP.inp);
	break;
      case WLZ_GREY_SHORT:
	gP.shp = baseGVP.shp + offset;
	*(gPP) = gP;
	(*(gVP)).shv = *(gP.shp);
	break;
      case WLZ_GREY_UBYTE:
	gP.ubp = baseGVP.ubp + offset;
	*(gPP) = gP;
	(*(gVP)).ubv = *(gP.ubp);
	break;
      case WLZ_GREY_FLOAT:
	gP.flp = baseGVP.flp + offset;
	*(gPP) = gP;
	(*(gVP)).flv = *(gP.flp);
	break;
      case WLZ_GREY_DOUBLE:
	gP.dbp = baseGVP.dbp + offset;
	*(gPP) = gP;
	(*(gVP)).dbv = *(gP.dbp);
	break;
      case WLZ_GREY_RGBA:
	gP.rgbp = baseGVP.rgbp + offset;
	*(gPP) = gP;
	(*(gVP)).rgbv = *(gP.rgbp);
	break;
      default:
	break;
    }
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Knowing that the given point is within the value table
*               computes the base pointer and offset for the point.
* \param	baseGVP			Destination pointer for the
*                                       base pointer.
* \param	offset			Destination pointer for the
*                                       offset from base pointer.
* \param	gVWSp			Grey value work space.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueComputeGreyP2D(WlzGreyP *baseGVP,
					   size_t *offset,
					   WlzGreyValueWSpace *gVWSp,
					   int line, int kol)
{
  int		itvCount,
  		kol0;
  WlzValueLine	*vLn;
  WlzValueIntervalLine *vILn;

  switch(gVWSp->gTabType2D)
  {
    case WLZ_GREY_TAB_RAGR:
      vLn = gVWSp->values2D.v->vtblines +
	    (line - gVWSp->values2D.v->line1);
      *baseGVP = vLn->values;
      *offset = kol - gVWSp->values2D.v->kol1 - vLn->vkol1;
      break;
    case WLZ_GREY_TAB_RECT:
      *baseGVP = gVWSp->values2D.r->values;
      *offset = (gVWSp->values2D.r->width *
                 (line - gVWSp->values2D.r->line1)) +
	        kol - gVWSp->values2D.r->kol1;
      break;
    case WLZ_GREY_TAB_INTL:
      kol0 = kol - gVWSp->values2D.i->kol1;
      vILn = gVWSp->values2D.i->vil + line - gVWSp->values2D.i->line1;
      itvCount = vILn->nintvs;
      vLn = vILn->vtbint;
      while(itvCount-- > 0)
      {
	if((kol0 >= vLn->vkol1 - 1) && (kol0 <= vLn->vlastkl))
	{
	  *baseGVP = vLn->values;
	  *offset = kol0 - vLn->vkol1;
	  itvCount = 0;				      /* Break from the loop */
	}
	++vLn;
      }
      break;
    case WLZ_GREY_TAB_TILED:
      WlzGreyValueComputeGreyPTiled2D(baseGVP, offset, gVWSp->values.t,
                                      line, kol);
      break;
    default:
      break;
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Knowing that the given point is within the value table
*               computes the base pointer and offset for the point
*               for a 2D tiled value table.
* \param	baseGVP			Destination pointer for the
*                                       base pointer.
* \param	offset			Destination pointer for the
*                                       offset from base pointer.
* \param	gVWSp			Grey value work space.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueComputeGreyPTiled2D(WlzGreyP *baseGVP,
				size_t *offset, WlzTiledValues *tVal,
				int line, int kol)
{
  WlzIVertex2 	rPos,
		tIdx;

  *offset = 0;
  (*baseGVP).v = NULL;
  rPos.vtX = kol - tVal->kol1;
  rPos.vtY = line - tVal->line1;
  tIdx.vtX = rPos.vtX / tVal->tileWidth;
  tIdx.vtY = rPos.vtY / tVal->tileWidth;
  if((tIdx.vtX >= 0) && (tIdx.vtX < tVal->nIdx[0]) &&
     (tIdx.vtY >= 0) && (tIdx.vtY < tVal->nIdx[1]))
  {
    size_t      idx;

    idx = (tIdx.vtY * tVal->nIdx[0]) + tIdx.vtX;
    idx = *(tVal->indices + idx);
    if((idx >= 0) && (idx < tVal->numTiles))
    {
      size_t    off;
      WlzIVertex2 tOff;

      tOff.vtX = rPos.vtX % tVal->tileWidth;
      tOff.vtY = rPos.vtY % tVal->tileWidth;
      off = (tOff.vtY * tVal->tileWidth) + tOff.vtX;
      (*baseGVP).v = tVal->tiles.v;
      *offset = (idx * tVal->tileSz) + off;
    }
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Knowing that the given point is within the value table
*               computes the base pointer and offset for the point
*               for a 3D tiled value table.
* \param	baseGVP			Destination pointer for the
*                                       base pointer.
* \param	offset			Destination pointer for the
*                                       offset from base pointer.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane coordinate of point.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueComputeGreyPTiled3D(WlzGreyP *baseGVP,
				size_t *offset, WlzTiledValues *tVal,
				int plane, int line, int kol)
{
  WlzIVertex3 	rPos,
		tIdx;

  *offset = 0;
  (*baseGVP).v = NULL;
  rPos.vtX = kol - tVal->kol1;
  rPos.vtY = line - tVal->line1;
  rPos.vtZ = plane - tVal->plane1;
  tIdx.vtX = rPos.vtX / tVal->tileWidth;
  tIdx.vtY = rPos.vtY / tVal->tileWidth;
  tIdx.vtZ = rPos.vtZ / tVal->tileWidth;
  if((tIdx.vtX >= 0) && (tIdx.vtX < tVal->nIdx[0]) &&
     (tIdx.vtY >= 0) && (tIdx.vtY < tVal->nIdx[1]) &&
     (tIdx.vtZ >= 0) && (tIdx.vtZ < tVal->nIdx[2]))
  {
    size_t	idx;

    idx = ((tIdx.vtZ * tVal->nIdx[1] + tIdx.vtY) * tVal->nIdx[0]) + tIdx.vtX;
    idx = *(tVal->indices + idx);
    if((idx >= 0) && (idx < tVal->numTiles))
    {
      size_t	off;
      WlzIVertex3 tOff;

      tOff.vtX = rPos.vtX % tVal->tileWidth;
      tOff.vtY = rPos.vtY % tVal->tileWidth;
      tOff.vtZ = rPos.vtZ % tVal->tileWidth;
      off = ((tOff.vtZ * tVal->tileWidth + tOff.vtY) * tVal->tileWidth) +
            tOff.vtX;
      (*baseGVP).v = tVal->tiles.v;
      *offset = (idx * tVal->tileSz) + off;
    }
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets a single grey value/pointer for the given point
*               from the 2D values and domain in the work space.
* \param	gVWSp			Grey value work space.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueGet2D1(WlzGreyValueWSpace *gVWSp,
				   int line, int kol)
{
  int		kol1,
		kolRel,
		count,
		valSet = 0;
  size_t 	offset;
  WlzGreyP	baseGVP;
  WlzInterval	*itv;
  WlzIntervalLine *itvLn;

  gVWSp->bkdFlag = 0;
  kol1 = gVWSp->iDom2D->kol1;
  if((gVWSp->values2D.core) && 
     (kol >= kol1) && (kol <= gVWSp->iDom2D->lastkl))
  {
    if((line >= gVWSp->iDom2D->line1) && (line <= gVWSp->iDom2D->lastln))
    {
      if(gVWSp->iDom2D->type == WLZ_INTERVALDOMAIN_RECT)
      {
	valSet = 1;
	WlzGreyValueComputeGreyP2D(&baseGVP, &offset, gVWSp, line, kol);
	WlzGreyValueSetGreyP(gVWSp->gVal, gVWSp->gPtr, gVWSp->gType,
			     baseGVP, offset);
      }
      else	          /* gVWSp->iDom2D->type == WLZ_INTERVALDOMAIN_INTVL */
      {
	itvLn = gVWSp->iDom2D->intvlines + line - gVWSp->iDom2D->line1;
	count = itvLn->nintvs;
	itv = itvLn->intvs;
	kolRel = kol - kol1;
	while(count-- > 0)
	{
	  if(kolRel < itv->ileft)
	  {
	    break;
	  }
	  else if(kolRel <= itv->iright)
	  {
	    valSet = 1;
	    WlzGreyValueComputeGreyP2D(&baseGVP, &offset, gVWSp, line, kol);
	    WlzGreyValueSetGreyP(gVWSp->gVal, gVWSp->gPtr, gVWSp->gType,
	    			 baseGVP, offset);
	    break;
	  }
	  ++itv;
	}
      }
    }
  }
  if(valSet == 0)
  {
    WlzGreyValueSetBkdP(gVWSp->gVal, gVWSp->gPtr, gVWSp->gType, gVWSp->gBkd);
    gVWSp->bkdFlag = 1;
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets a single grey value/pointer for the given point
*               from the tiled 3D values and domain in the work space.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane coordinate of point.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueGet3DTiled(WlzGreyValueWSpace *gVWSp,
				       int plane, int line, int kol)
{
  int           valSet = 0;

  int           plRel;
  WlzDomain     *dom;

  gVWSp->bkdFlag = 0;
  plRel = plane - gVWSp->domain.p->plane1;
  if((plRel >= 0) && (plane <= gVWSp->domain.p->lastpl))
  {
    dom = gVWSp->domain.p->domains + plRel;
    if((dom != NULL) && ((*dom).core != NULL))
    {
      int	kol1;
      WlzIntervalDomain *iDom;

      iDom = (*dom).i;
      kol1 = iDom->kol1;
      if((kol >= kol1) && (kol <= iDom->lastkl))
      {
	if((line >= iDom->line1) && (line <= iDom->lastln))
	{
	  if(iDom->type == WLZ_INTERVALDOMAIN_RECT)
	  {
	    valSet = 1;
	  }
	  else
	  {
	    int count,
		kolRel,
		lnRel;
	    WlzInterval *itv;
	    WlzIntervalLine *itvLn;

	    lnRel = line - iDom->line1;
	    itvLn = iDom->intvlines + lnRel;
	    count = itvLn->nintvs;
	    itv = itvLn->intvs;
	    kolRel = kol - kol1;
	    while(count-- > 0)
	    {
	      if(kolRel < itv->ileft)
	      {
		break;
	      }
	      else if(kolRel <= itv->iright)
	      {
		valSet = 1;
		break;
	      }
	      ++itv;
	    }
	  }
	  if(valSet)
	  {
	    size_t   	offset;
	    WlzGreyP 	baseGVP;

	    WlzGreyValueComputeGreyPTiled3D(&baseGVP, &offset,
					    gVWSp->values.t,
					    plane, line, kol);
	    WlzGreyValueSetGreyP(gVWSp->gVal, gVWSp->gPtr, gVWSp->gType,
				 baseGVP, offset);
	  }
	}
      }
    }
  }
  if(valSet == 0)
  {
    WlzGreyValueSetBkdP(gVWSp->gVal, gVWSp->gPtr, gVWSp->gType, gVWSp->gBkd);
    gVWSp->bkdFlag = 1;
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets a single grey value/pointer for the given point
*               from the 3D values and domain in the work space.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane coordinate of point.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueGet3D1(WlzGreyValueWSpace *gVWSp,
				   int plane, int line, int kol)
{
  int           plRel,
                valSet = 0;
  WlzDomain     *domP;

  gVWSp->bkdFlag = 0;
  plRel = plane - gVWSp->domain.p->plane1;
  if((plRel >= 0) && (plane <= gVWSp->domain.p->lastpl))
  {
    if(gVWSp->plane == plane)
    {
      WlzGreyValueGet2D1(gVWSp, line, kol);
      valSet = 1;
    }
    else
    {
      WlzValues *valP;

      domP = gVWSp->domain.p->domains + plRel;
      valP = gVWSp->values.vox->values + plRel;
      /* Check for non-NULL domain and valuetable pointers until empty
       * obj consistently implemented. */
      if(domP && valP && (*domP).core && (*valP).core)
      {
	gVWSp->plane = plane;
	gVWSp->iDom2D = (*domP).i;
	gVWSp->values2D = (*valP);
	gVWSp->gTabType2D = gVWSp->gTabTypes3D[plRel];
	WlzGreyValueGet2D1(gVWSp, line, kol);
	valSet = 1;
      }
    }
  }
  if(valSet == 0)
  {
    WlzGreyValueSetBkdP(gVWSp->gVal, gVWSp->gPtr, gVWSp->gType, gVWSp->gBkd);
    gVWSp->bkdFlag = 1;
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets four grey values/pointers for the given point
*               from the 2D values and domain in the work space.
* \param	gVWSp			Grey value work space.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueGet2DCon(WlzGreyValueWSpace *gVWSp,
				     int line, int kol)
{
  int		kol1,
		kolRel,
		count,
  		pass;
  size_t 	offset;
  unsigned	hitMsk;
  WlzGreyP	baseGVP;
  WlzInterval	*itv;
  WlzIntervalLine *itvLn;
  WlzGreyV	*gVP;
  WlzGreyP	*gPP;

  hitMsk = 0;
  gVWSp->bkdFlag = 0;
  kol1 = gVWSp->iDom2D->kol1;
  gVP = gVWSp->gVal;
  gPP = gVWSp->gPtr;
  WlzGreyValueSetBkdPN(gVP, gPP, gVWSp->gType, gVWSp->gBkd, 4);
  if((gVWSp->values2D.core) &&
     ((kol + 1) >= kol1) && (kol <= gVWSp->iDom2D->lastkl))
  {
    pass = 0;
    while(pass < 2)
    {
      if((line >= gVWSp->iDom2D->line1) && (line <= gVWSp->iDom2D->lastln))
      {
	if(gVWSp->iDom2D->type == WLZ_INTERVALDOMAIN_RECT)
	{
	  WlzGreyValueComputeGreyP2D(&baseGVP, &offset, gVWSp, line, kol);
	  if(kol >= kol1)
	  {
	    hitMsk |= 1 << (pass * 2);
	    WlzGreyValueSetGreyP(gVP, gPP, gVWSp->gType, baseGVP, offset);
	  }
	  if(kol <= (gVWSp->iDom2D->lastkl - 1))
	  {
	    hitMsk |= 1 << ((pass * 2) + 1);
	    WlzGreyValueSetGreyP(gVP + 1, gPP + 1, gVWSp->gType,
	    			 baseGVP, offset + 1);
	  }
	}
	else		  /* gVWSp->iDom2D->type == WLZ_INTERVALDOMAIN_INTVL */
	{
	  itvLn = gVWSp->iDom2D->intvlines + line - gVWSp->iDom2D->line1;
	  count = itvLn->nintvs;
	  itv = itvLn->intvs;
	  kolRel = kol - kol1;
	  while(count-- > 0)
	  {
	    if((kolRel + 1) < itv->ileft)
	    {
	      count = 0;			      /* Break from the loop */
	    }
	    else if(kolRel <= itv->iright)
	    {
	      WlzGreyValueComputeGreyP2D(&baseGVP, &offset, gVWSp, line, kol);
	      if(kol >= (kol1 + itv->ileft))
	      {
	        hitMsk |= 1 << (pass * 2);
		WlzGreyValueSetGreyP(gVP, gPP, gVWSp->gType, baseGVP, offset);
	      }
	      if(kol <= (itv->iright + kol1 - 1))
	      {
	        hitMsk |= 1 << ((pass * 2) + 1);
		WlzGreyValueSetGreyP(gVP + 1, gPP + 1, gVWSp->gType,
				     baseGVP, offset + 1);
	      }
	    }
	    ++itv;
	  }
	}
      }
      ++pass;
      ++line;
      gVP += 2;
      gPP += 2;
    }
  }
  gVWSp->bkdFlag = 0xf & ~hitMsk;
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets eight grey value/pointers for the given point
*               from the 3D values and domain in the work space.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane coordinate of point.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueGet3DCon(WlzGreyValueWSpace *gVWSp,
				   int plane, int line, int kol)
{
  int		tI0,
  		planeOff,
		planeRel,
		savePlane;
  WlzDomain	*domP;
  WlzValues	*valP;
  WlzObjectType	saveGTabType2D;
  WlzIntervalDomain *saveIDom2D;
  WlzValues	saveValues2D;
  int		planeSet[2];
  WlzGreyP	saveGPtr[4];
  WlzGreyV	saveGVal[4];

  planeOff = 0;
  while(planeOff < 2)
  {
    planeSet[planeOff] = 0;
    if((plane > gVWSp->domain.p->plane1) && (plane < gVWSp->domain.p->lastpl))
    {
      if(plane == gVWSp->plane)
      {
	WlzGreyValueGet2DCon(gVWSp, line, kol);
	planeSet[planeOff] = 1;
      }
      else
      {
	planeRel = plane - gVWSp->domain.p->plane1;
	domP = gVWSp->domain.p->domains + planeRel;
	valP = gVWSp->values.vox->values + planeRel;
	if(domP && valP)
	{
          if(planeSet[0])
	  {
	    savePlane = gVWSp->plane;
	    saveIDom2D = gVWSp->iDom2D;
	    saveValues2D = gVWSp->values2D;
	    saveGTabType2D = gVWSp->gTabType2D;
	  }
	  gVWSp->plane = plane;
	  gVWSp->iDom2D = (*domP).i;
	  gVWSp->values2D = (*valP);
	  gVWSp->gTabType2D = gVWSp->gTabTypes3D[planeRel];
	  WlzGreyValueGet2D1(gVWSp, line, kol);
	  planeSet[planeOff] = 1;
	}
      }
    }
    if(planeOff == 0)
    {
      if(planeSet[0])
      {
	tI0 = 4;
	while(--tI0 >= 0)
	{
	  saveGPtr[tI0] = gVWSp->gPtr[tI0];
	  saveGVal[tI0] = gVWSp->gVal[tI0];
	}
      }
      else
      {
        WlzGreyValueSetBkdPN(saveGVal, saveGPtr,
			     gVWSp->gType, gVWSp->gBkd, 4);
      }
    }
    ++plane;
    ++planeOff;
  }
  if(planeSet[0] && planeSet[1])
  {
    gVWSp->plane = savePlane;
    gVWSp->iDom2D = saveIDom2D;
    gVWSp->values2D = saveValues2D;
    gVWSp->gTabType2D = saveGTabType2D;
  }
  tI0 = 4;
  while(--tI0 >= 0)
  {
    gVWSp->gPtr[tI0 + 4] = gVWSp->gPtr[tI0];
    gVWSp->gPtr[tI0] = saveGPtr[tI0];
    gVWSp->gVal[tI0 + 4] = gVWSp->gVal[tI0];
    gVWSp->gVal[tI0] = saveGVal[tI0];
  }
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets eight grey value/pointers for the given point
*               from the 3D tiled values and domain in the work space.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane coordinate of point.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueGet3DConTiled(WlzGreyValueWSpace *gVWSp,
				   int plane, int line, int kol)
{
  int		idK,
  		idL,
		idP;
  unsigned	valMsk = 0;

  /* First set bit mask for which voxels are inside the domain. */
  for(idP = 0; idP < 2; ++idP)
  {
    int		idV,
    		pl,
    		plRel;

    pl = plane + idP;
    plRel = pl - gVWSp->domain.p->plane1;
    if((plRel >= 0) && (pl <= gVWSp->domain.p->lastpl))
    {
      WlzDomain *dom;

      dom = gVWSp->domain.p->domains + plRel;
      if((dom != NULL) && ((*dom).i != NULL))
      {
        WlzIntervalDomain *iDom;

        iDom = (*dom).i;
	for(idL = 0; idL < 2; ++idL)
	{
	  int	ln;

	  ln = line + idL;
	  if((ln >= iDom->line1) && (ln <= iDom->line1))
	  {
	    if((kol + 1 >= iDom->kol1) &&
	       (kol <= iDom->lastkl))
	    {
	      idV = ((idP << 1) + idL) << 1;
	      if(iDom->type == WLZ_INTERVALDOMAIN_RECT)
	      {
	        valMsk |= ((kol >= iDom->kol1) |
		          ((kol < iDom->lastkl) << 1)) << idV;

	      }
	      else /* iDom->type == WLZ_INTERVALDOMAIN_INTVL */
	      {
		int	idI,
			klRel,
			lnRel;
		WlzIntervalLine *itvLn;
		WlzInterval *itv;

	        lnRel = ln - iDom->line1;
		itvLn = iDom->intvlines + lnRel;
		itv = itvLn->intvs;
		klRel = kol - iDom->kol1;
		for(idI = 0; idI < itvLn->nintvs; ++idI)
		{
		  if(klRel + 1 < itv->ileft)
		  {
		    break;
		  }
		  else if(klRel <= itv->iright)
		  {
		    valMsk |= ((klRel >= itv->ileft) |
		               ((klRel < itv->iright) << 1)) << idV;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  /* Now set grey values and pointers knowing which are within the domain. */
  if(gVWSp->bkdFlag == 0xff)
  {
    WlzGreyValueSetBkdPN(gVWSp->gVal, gVWSp->gPtr,
                         gVWSp->gType, gVWSp->gBkd, 8);
  }
  else
  {
    int		idV;
    WlzIVertex3	rPos,
    		tIdx,
		tOff;
    size_t 	offset;
    WlzTiledValues *tVal;

    idV = 0;
    tVal = gVWSp->values.t;
    for(idP = 0; idP < 2; ++idP)
    {
      rPos.vtZ = plane - tVal->plane1 + idP;
      tIdx.vtZ = (rPos.vtZ / tVal->tileWidth) * tVal->nIdx[1];
      tOff.vtZ = (rPos.vtZ % tVal->tileWidth) * tVal->tileWidth;
      for(idL = 0; idL < 2; ++idL)
      {
        rPos.vtY = line - tVal->line1 + idL;
	tIdx.vtY = (tIdx.vtZ + (rPos.vtY / tVal->tileWidth)) * tVal->nIdx[0];
	tOff.vtY = (tOff.vtZ + (rPos.vtY % tVal->tileWidth)) * tVal->tileWidth;
	for(idK = 0; idK < 2; ++idK)
	{
	  if((valMsk & (1 < idV)) == 0)
	  {
	    WlzGreyValueSetBkdPN(gVWSp->gVal + idV, gVWSp->gPtr + idV,
	                         gVWSp->gType, gVWSp->gBkd, 1);
	  }
	  else
	  {
            rPos.vtX = line - tVal->kol1 + idK;
	    tIdx.vtX = tIdx.vtY + (rPos.vtX / tVal->tileWidth);
            tOff.vtX = tOff.vtY + (rPos.vtX % tVal->tileWidth);
            offset = *(tVal->indices + tIdx.vtX) * tVal->tileSz + tOff.vtX;
	    WlzGreyValueSetGreyP(gVWSp->gVal + idV, gVWSp->gPtr + idV,
	                         gVWSp->gType, tVal->tiles, offset);
	  }
	  ++idV;
	}
      }
    }
  }
  gVWSp->bkdFlag = 0xff & (~valMsk);
}

/*!
* \return	void
* \ingroup	WlzAccess
* \brief	Gets four or eight grey value/pointers for the given
*               point from the 2D or 3D values and domain in the work
*               space.
* \param	gVWSp			Grey value work space.
* \param	plane			Plane coordinate of point.
* \param	line			Line coordinate of point.
* \param	kol			Column coordinate of point.
*/
static void	WlzGreyValueGetTransCon(WlzGreyValueWSpace *gVWSp,
				        int plane, int line, int kol)
{
  int		idN,
  		idX,
  		idY,
		idZ;
  WlzDVertex2	vtx2;
  WlzDVertex3	vtx3;
  WlzGreyP	gPtr[8];
  WlzGreyV	gVal[8];

  switch(gVWSp->objType)
  {
    case WLZ_2D_DOMAINOBJ:
      idN = 0;
      for(idY = 0; idY < 2; ++idY)
      {
	for(idX = 0; idX < 2; ++idX)
	{
	  vtx2.vtX = kol + idX;
	  vtx2.vtY = line + idY;
	  vtx2 = WlzAffineTransformVertexD2(gVWSp->invTrans, vtx2, NULL);
	  WlzGreyValueGet2D1(gVWSp, (int )(vtx2.vtY), (int )(vtx2.vtX));
	  gPtr[idN] = gVWSp->gPtr[0];
	  gVal[idN] = gVWSp->gVal[0];
	  ++idN;
	}
      }
      for(idN = 0; idN < 4; ++idN)
      {
        gVWSp->gPtr[idN] = gPtr[idN];
	gVWSp->gVal[idN] = gVal[idN];
      }
      break;
    case WLZ_3D_DOMAINOBJ:
      idN = 0;
      for(idZ = 0; idZ < 2; ++idZ)
      {
	for(idY = 0; idY < 2; ++idY)
	{
	  for(idX = 0; idX < 2; ++idX)
	  {
	    vtx3.vtX = kol + idX;
	    vtx3.vtY = line + idY;
	    vtx3.vtZ = plane + idZ;
	    vtx3 = WlzAffineTransformVertexD3(gVWSp->invTrans, vtx3, NULL);
	    if(gVWSp->gTabType == WLZ_GREY_TAB_TILED)
	    {
	      WlzGreyValueGet3DTiled(gVWSp,
			 (int )(vtx3.vtZ), (int )(vtx3.vtY), (int )(vtx3.vtX));
	    }
	    else
	    {
	      WlzGreyValueGet3D1(gVWSp,
			 (int )(vtx3.vtZ), (int )(vtx3.vtY), (int )(vtx3.vtX));
	    }
	    gPtr[idN] = gVWSp->gPtr[0];
	    gVal[idN] = gVWSp->gVal[0];
	    ++idN;
	  }
	}
      }
      for(idN = 0; idN < 8; ++idN)
      {
        gVWSp->gPtr[idN] = gPtr[idN];
	gVWSp->gVal[idN] = gVal[idN];
      }
      break;
    default:
      break;
  }
}
