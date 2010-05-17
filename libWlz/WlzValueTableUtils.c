#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzValueTableUtils_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzValueTableUtils.c
* \author       Bill Hill, Richard Baldock
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
* \brief	Functions for computing value amd value table types.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Type of grey table.
* \ingroup	WlzValuesUtils
* \brief	Computes a grey table type from table and grey types.
* \param	tableType		The basic table type.
* \param	greyType		The grey type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType	WlzGreyTableType(WlzObjectType tableType,
			         WlzGreyType greyType,
				 WlzErrorNum *dstErr)
{
  WlzObjectType gTabType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(tableType)
  {
    case WLZ_GREY_TAB_RAGR:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_RAGR_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_RAGR_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_RAGR_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_RAGR_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_RAGR_DOUBLE;
	  break;
        case WLZ_GREY_RGBA:
	  gTabType = WLZ_VALUETABLE_RAGR_RGBA;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_RECT:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_RECT_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_RECT_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_RECT_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_RECT_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_RECT_DOUBLE;
	  break;
        case WLZ_GREY_RGBA:
	  gTabType = WLZ_VALUETABLE_RECT_RGBA;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_INTL:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_INTL_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_INTL_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_INTL_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_INTL_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_INTL_DOUBLE;
	  break;
        case WLZ_GREY_RGBA:
	  gTabType = WLZ_VALUETABLE_INTL_RGBA;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_TILED:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_TILED_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_TILED_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_TILED_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_TILED_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_TILED_DOUBLE;
	  break;
        case WLZ_GREY_RGBA:
	  gTabType = WLZ_VALUETABLE_TILED_RGBA;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    default:
      errNum = WLZ_ERR_VALUES_TYPE;
      break;
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(gTabType);
}

/*!
* \return	Type of grey value.
* \ingroup	WlzValuesUtils
* \brief	Computes the type of grey from a grey table type.
* \param	gTabType		The basic table type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGreyType	WlzGreyTableTypeToGreyType(WlzObjectType gTabType,
				           WlzErrorNum *dstErr)
{
  WlzGreyType	greyType = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gTabType)
  {
    case WLZ_VALUETABLE_RAGR_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_INT:
      greyType = WLZ_GREY_INT;
      break;
    case WLZ_VALUETABLE_RAGR_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_SHORT:
      greyType = WLZ_GREY_SHORT;
      break;
    case WLZ_VALUETABLE_RAGR_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_UBYTE:
      greyType = WLZ_GREY_UBYTE;
      break;
    case WLZ_VALUETABLE_RAGR_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_FLOAT:
      greyType = WLZ_GREY_FLOAT;
      break;
    case WLZ_VALUETABLE_RAGR_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_DOUBLE:
      greyType = WLZ_GREY_DOUBLE;
      break;
    case WLZ_VALUETABLE_RAGR_RGBA: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_RGBA: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_RGBA: /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_RGBA:
      greyType = WLZ_GREY_RGBA;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(greyType);
}

/*!
* \return	Type of grey value.
* \ingroup	WlzValuesUtils
* \brief	Computes the type of table from a  grey table type.
* \param	gTabType		The basic table type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType WlzGreyTableTypeToTableType(WlzObjectType gTabType,
				 	  WlzErrorNum *dstErr)
{
  WlzObjectType	tableType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gTabType)
  {
    case WLZ_VALUETABLE_RAGR_INT:    /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_SHORT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_UBYTE:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_FLOAT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_RGBA:
      tableType = WLZ_GREY_TAB_RAGR;
      break;
    case WLZ_VALUETABLE_RECT_INT:    /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_SHORT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_UBYTE:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_FLOAT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_DOUBLE  /* FALLTHROUGH */:
    case WLZ_VALUETABLE_RECT_RGBA:
      tableType = WLZ_GREY_TAB_RECT;
      break;
    case WLZ_VALUETABLE_INTL_INT:    /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_SHORT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_UBYTE:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_FLOAT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_RGBA:
      tableType = WLZ_GREY_TAB_INTL;
      break;
    case WLZ_VALUETABLE_TILED_INT:    /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_SHORT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_UBYTE:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_FLOAT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_RGBA:
      tableType = WLZ_GREY_TAB_TILED;
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tableType);
}

/*!
* \return	WLZ_GREY_TAB_TILED if the grey table type
* 		is tiled otherwise zero.
* \ingroup	WlzValuesUtils
* \brief	Determines whether the grey table type is tiled.
* \param	gTabType		The basic table type.
*/
WlzObjectType	WlzGreyTableIsTiled(WlzObjectType gTabType)
{
  WlzObjectType	is = (WlzObjectType )0;

  switch(gTabType)
  {
    case WLZ_VALUETABLE_TILED_INT:    /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_SHORT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_UBYTE:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_FLOAT:  /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_TILED_RGBA:
      is = WLZ_GREY_TAB_TILED;
      break;
    default:
      break;
  }
  return(is);
}

/*!
* \return	Grey type, WLZ_GREY_ERROR on error.
* \ingroup	WlzValuesUtils
* \brief	Gets the grey type of the values in a Woolz object.
*               If the object is not a domain object with values
*               an error is returned. If the object is a 3D domain
*               object with values, all 2D value tables are checked
*               and an error is returned if they don't all have the
*               same grey type.
* \param	obj			Given object.
* \param	dstErr			Destination ptr for error code,
*					may be NULL.
*/
WlzGreyType	WlzGreyTypeFromObj(WlzObject *obj, WlzErrorNum *dstErr)
{
  int		pCnt,
    pIdx,
    firstDomFlg;
  WlzDomain	dom2D;
  WlzValues	val2D;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(obj->type)
    {
    case WLZ_2D_DOMAINOBJ:
      gType = WlzGreyTableTypeToGreyType(obj->values.core->type, &errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
      if(WlzGreyTableIsTiled(obj->values.core->type))
      {
        gType = WlzGreyTableTypeToGreyType(obj->values.core->type, &errNum);
      }
      else
      {
	pIdx = 0;
	firstDomFlg = 1;
	pCnt = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
	while((errNum == WLZ_ERR_NONE) && (pIdx < pCnt))
	{
	  if(((dom2D = *(obj->domain.p->domains + pIdx)).core != NULL) &&
	     (dom2D.core->type != WLZ_EMPTY_DOMAIN))
	  {
	    if(((val2D = *(obj->values.vox->values + pIdx)).core == NULL) ||
	       (val2D.core->type == WLZ_EMPTY_VALUES))
	    {
	      errNum = WLZ_ERR_VALUES_NULL;
	    }
	    else
	    {
	      if(firstDomFlg)
	      {
		gType = WlzGreyTableTypeToGreyType(val2D.core->type, &errNum);
	      }
	      else
	      {
		if(gType != WlzGreyTableTypeToGreyType(val2D.core->type,
						       &errNum))
		{
		  errNum = WLZ_ERR_VALUES_DATA;
		}
	      }
	    }
	  }
	  ++pIdx;
	}
      }
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(gType);
}

/*!
* \return	Voxel size.
* \ingroup	WlzValuesUtils
* \brief	Gets the given 3D domain objects voxel size.
* \param	obj			Given 3D domain object.
* \param	dstErr			Destination pointer for error code,
* 					may be NULL.
*/
WlzDVertex3	WlzVozelSz(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzDVertex3	sz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  sz.vtX = sz.vtY = sz.vtZ = 0;
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    sz.vtX = obj->domain.p->voxel_size[0];
    sz.vtY = obj->domain.p->voxel_size[1];
    sz.vtZ = obj->domain.p->voxel_size[2];
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(sz);
}

/*!
* \return	Value pointer which may be null on error.
* \ingroup	WlzValuesUtils
* \brief	Gets a pointer to the valuetable entry for the given
* 		index. If the indexed values are not valid for the
* 		given index NULL will be returned. See also
* 		WlzIndexedValueExtGet().
* \param	ixv			Indexed value table.
* \param	idx			Given index.
*/
void		*WlzIndexedValueGet(WlzIndexedValues *ixv, int idx)
{
  void 		*val = NULL;

  if((ixv != NULL) && (idx >= 0))
  {
    val = AlcVectorItemGet(ixv->values, idx);
  }
  return(val);
}

/*!
* \return	Value pointer which may be null on error.
* \ingroup	WlzValuesUtils
* \brief	Gets a pointer to the valuetable entry for the given
* 		index. The value table is extended as required so that
* 		there should always be a valid entry unless memory
* 		allocation fails. See also WlzIndexedValueGet().
* \param	ixv			Indexed value table.
* \param	idx			Given index.
*/
void		*WlzIndexedValueExtGet(WlzIndexedValues *ixv, int idx)
{
  void 		*val = NULL;

  if((ixv != NULL) && (idx >= 0))
  {
    val = AlcVectorExtendAndGet(ixv->values, idx);
  }
  return(val);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	Sets all values covered by the domain in the given object
* 		to the given value. Ths function calls memcpy() for each
* 		indexed value using the given value count and pointer.
* \param	obj			Given object which must be a CMesh
* 					object with indexed values.
* \param	cnt			Size of the values to be copied.
* \param	val			Basic value set for the indexed
* 					values. Type must be the same as
* 					the indexed value table value type.
*/
WlzErrorNum     WlzIndexedValuesSet(WlzObject *obj, size_t cnt, void *val)
{
  int		idx,
  		max;
  void		*dst;
  AlcVector	*vec;
  WlzCMeshElm2D	*e2;
  WlzCMeshElm2D5 *e2d5;
  WlzCMeshElm3D	*e3;
  WlzCMeshNod2D	*n2;
  WlzCMeshNod2D5 *n2d5;
  WlzCMeshNod3D	*n3;
  WlzCMeshP	mesh;
  WlzIndexedValues *ixv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((ixv = obj->values.x) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(ixv->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_CMESH_2D:
	mesh.m2 = obj->domain.cm2;
	if(mesh.m2->type != WLZ_CMESH_2D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  switch(ixv->attach)
	  {
	    case WLZ_VALUE_ATTACH_NOD:
	      max = mesh.m2->res.nod.maxEnt;
	      if(WlzIndexedValueExtGet(ixv, max) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
	        vec = mesh.m2->res.nod.vec;
		for(idx = 0; idx < max; ++idx)
		{
		  n2 = (WlzCMeshNod2D *)AlcVectorItemGet(vec, idx);
		  if(n2->idx >= 0)
		  {
		    dst = WlzIndexedValueGet(ixv, idx);
		    memcpy(dst, val, cnt);
		  }
		}
	      }
	      break;
	    case WLZ_VALUE_ATTACH_ELM:
	      max = mesh.m2->res.elm.maxEnt;
	      if(WlzIndexedValueExtGet(ixv, max) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
	        vec = mesh.m2->res.elm.vec;
		for(idx = 0; idx < max; ++idx)
		{
		  e2 = (WlzCMeshElm2D *)AlcVectorItemGet(vec, idx);
		  if(e2->idx >= 0)
		  {
		    dst = WlzIndexedValueGet(ixv, idx);
		    memcpy(dst, val, cnt);
		  }
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_VALUES_DATA;
	      break;
	  }
	}
        break;
      case WLZ_CMESH_2D5:
	mesh.m2d5 = obj->domain.cm2d5;
	if(mesh.m2->type != WLZ_CMESH_2D5)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  switch(ixv->attach)
	  {
	    case WLZ_VALUE_ATTACH_NOD:
	      max = mesh.m2d5->res.nod.maxEnt;
	      if(WlzIndexedValueExtGet(ixv, max) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
	        vec = mesh.m2d5->res.nod.vec;
		for(idx = 0; idx < max; ++idx)
		{
		  n2d5 = (WlzCMeshNod2D5 *)AlcVectorItemGet(vec, idx);
		  if(n2d5->idx >= 0)
		  {
		    dst = WlzIndexedValueGet(ixv, idx);
		    memcpy(dst, val, cnt);
		  }
		}
	      }
	      break;
	    case WLZ_VALUE_ATTACH_ELM:
	      max = mesh.m2d5->res.elm.maxEnt;
	      if(WlzIndexedValueExtGet(ixv, max) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
	        vec = mesh.m2d5->res.elm.vec;
		for(idx = 0; idx < max; ++idx)
		{
		  e2d5 = (WlzCMeshElm2D5 *)AlcVectorItemGet(vec, idx);
		  if(e2d5->idx >= 0)
		  {
		    dst = WlzIndexedValueGet(ixv, idx);
		    memcpy(dst, val, cnt);
		  }
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_VALUES_DATA;
	      break;
	  }
	}
        break;
      case WLZ_CMESH_3D:
	mesh.m3 = obj->domain.cm3;
	if(mesh.m3->type != WLZ_CMESH_3D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  switch(ixv->attach)
	  {
	    case WLZ_VALUE_ATTACH_NOD:
	      max = mesh.m3->res.nod.maxEnt;
	      if(WlzIndexedValueExtGet(ixv, max) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
	        vec = mesh.m3->res.nod.vec;
		for(idx = 0; idx < max; ++idx)
		{
		  n3 = (WlzCMeshNod3D *)AlcVectorItemGet(vec, idx);
		  if(n3->idx >= 0)
		  {
		    dst = WlzIndexedValueGet(ixv, idx);
		    memcpy(dst, val, cnt);
		  }
		}
	      }
	      break;
	    case WLZ_VALUE_ATTACH_ELM:
	      max = mesh.m3->res.elm.maxEnt;
	      if(WlzIndexedValueExtGet(ixv, max) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
	        vec = mesh.m3->res.elm.vec;
		for(idx = 0; idx < max; ++idx)
		{
		  e3 = (WlzCMeshElm3D *)AlcVectorItemGet(vec, idx);
		  if(e3->idx >= 0)
		  {
		    dst = WlzIndexedValueGet(ixv, idx);
		    memcpy(dst, val, cnt);
		  }
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_VALUES_DATA;
	      break;
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  return(errNum);
}

