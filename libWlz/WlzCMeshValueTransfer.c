#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshValueTransfer_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzCMeshValueTransfer.c
* \author       Bill Hill
* \date         August 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Functions for transfering values from one mesh to
* 		another.
* \ingroup	WlzMesh
*/

#include <Wlz.h>

static WlzObject		*WlzCMeshValueTransferD2D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshValueTransferM2D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshValueTransferD3D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshValueTransferM3D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);

/*!
* \return	New mesh object or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the mesh domain of the target
* 		object and (within the domain intersection) the values
* 		of the source object. Given a pair of 2, 2.5 or 3D
* 		conforming meshes, with both of the same type, a new
* 		mesh object is created using the domain of the target
*               and with the value type of the source object. Values are
*               transfered from the source to the target mesh within their
*               intersection. Outside of the intersection values are set
*               to the given external value.
* \param	srcObj			Given source mesh object which must
* 					be of the same type as the target and
* 					have values attached.
* \param	tgtObj			Given target mesh object which must
* 					be of the same type as the source.
* \param	extVal			Given external value to use outside
* 					of the source - target mesh
* 					intersection.
* \param	interp			Interpolation method.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzCMeshValueTransfer(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((srcObj == NULL) || (tgtObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((srcObj->domain.core == NULL) || (tgtObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(interp)
    {
      case WLZ_INTERPOLATION_NEAREST: /* FALLTHROUGH */
      case WLZ_INTERPOLATION_LINEAR:
        break;
      case WLZ_INTERPOLATION_BARYCENTRIC:
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      default:
        errNum = WLZ_ERR_PARAM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(tgtObj->type)
    {
      case WLZ_CMESH_2D:
	switch(srcObj->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	    rtnObj = WlzCMeshValueTransferD2D(srcObj, tgtObj, extVal, interp,
					       &errNum);
	    break;
	  case WLZ_CMESH_2D:
	    rtnObj = WlzCMeshValueTransferM2D(srcObj, tgtObj, extVal, interp,
					       &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
        break;
      case WLZ_CMESH_2D5:
        errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	switch(srcObj->type)
	{
	  case WLZ_3D_DOMAINOBJ:
	    rtnObj = WlzCMeshValueTransferD3D(srcObj, tgtObj, extVal, interp,
					       &errNum);
	    break;
	  case WLZ_CMESH_3D:
	    rtnObj = WlzCMeshValueTransferM3D(srcObj, tgtObj, extVal, interp,
					       &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
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
  return(rtnObj);
}

/*!
* \return	New mesh object or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the mesh domain of the target
* 		object and (within the domain intersection) the values
* 		of the source spatial domain object. Given a pair of 2D
* 		spatial domain and conforming mesh objects. The resulting
* 		mesh values will always be WLZ_GREY_DOUBLE and attached
* 		to the mesh nodes.  See WlzCMeshValueTransfer().
* \param	srcObj			Given source spatial domain object.
* \param	tgtObj			Given target mesh object.
* \param	extVal			Given external value.
* \param	interp			Interpolation method.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzCMeshValueTransferD2D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr)
{
  WlzCMesh2D	   *rtnMesh;
  WlzIndexedValues *rtnIxv = NULL;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gVWSp->gType)
    {
      case WLZ_GREY_INT:    /* FALLTHROUGH */
      case WLZ_GREY_SHORT:  /* FALLTHROUGH */
      case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
      case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
        switch(interp)
	{
	  case WLZ_INTERPOLATION_NEAREST: /* FALLTHROUGH */
	  case WLZ_INTERPOLATION_LINEAR:
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_TYPE;
	    break;
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rtnMesh = tgtObj->domain.cm2;
    errNum = WlzValueConvertPixel(&extVal, extVal, WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rtnIxv = WlzMakeIndexedValues(tgtObj, 0, NULL,
				  WLZ_GREY_DOUBLE, WLZ_VALUE_ATTACH_NOD,
				  &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	rtnVal;

    rtnVal.x = rtnIxv;
    rtnObj = WlzMakeMain(tgtObj->type, tgtObj->domain, rtnVal, NULL, NULL,
                         &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		   rNI;
    WlzCMeshEntRes *rNR;

    rNR = &(rtnMesh->res.nod);
    for(rNI = 0; rNI < rNR->maxEnt; ++rNI)
    {
      WlzCMeshNod2D *rNod;

      rNod = (WlzCMeshNod2D *)AlcVectorItemGet(rNR->vec, rNI);
      if(rNod->idx >= 0)
      {
	WlzGreyP rGP;

	rGP.v = WlzIndexedValueGet(rtnIxv, rNI);
	switch(interp)
	{
	  case WLZ_INTERPOLATION_NEAREST:
	    WlzGreyValueGet(gVWSp, 0, rNod->pos.vtY, rNod->pos.vtX);
	    if(gVWSp->bkdFlag)
	    {
	      rGP.dbp[0] = extVal.v.dbv;
	    }
	    else
	    {
	      switch(gVWSp->gType)
	      {
	        case WLZ_GREY_INT:
		  rGP.dbp[0] = (gVWSp->gVal[0]).inv;
		  break;
	        case WLZ_GREY_SHORT:
		  rGP.dbp[0] = (gVWSp->gVal[0]).shv;
		  break;
	        case WLZ_GREY_UBYTE:
		  rGP.dbp[0] = (gVWSp->gVal[0]).ubv;
		  break;
	        case WLZ_GREY_FLOAT:
		  rGP.dbp[0] = (gVWSp->gVal[0]).flv;
		  break;
	        case WLZ_GREY_DOUBLE:
		  rGP.dbp[0] = (gVWSp->gVal[0]).dbv;
		  break;
	        default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	    }
	    break;
	  case WLZ_INTERPOLATION_LINEAR:
	    WlzGreyValueGetCon(gVWSp, 0, rNod->pos.vtY, rNod->pos.vtX);
	    if(gVWSp->bkdFlag)          
	    {
	      rGP.dbp[0] = extVal.v.dbv;
	    }
	    else
	    {
	      int	i;
	      double	iBuf[4],
	      		vBuf[4];

	      switch(gVWSp->gType)
	      {
	        case WLZ_GREY_INT:
		  for(i = 0; i < 4; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).inv;
		  }
		  break;
	        case WLZ_GREY_SHORT:
		  for(i = 0; i < 4; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).shv;
		  }
		  break;
	        case WLZ_GREY_UBYTE:
		  for(i = 0; i < 4; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).ubv;
		  }
		  break;
	        case WLZ_GREY_FLOAT:
		  for(i = 0; i < 4; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).flv;
		  }
		  break;
	        case WLZ_GREY_DOUBLE:
		  for(i = 0; i < 4; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).dbv;
		  }
		  break;
	        default:
		  break;
	      }
	      iBuf[0] = rNod->pos.vtX - WLZ_NINT(rNod->pos.vtX - 0.5);
	      iBuf[1] = rNod->pos.vtY - WLZ_NINT(rNod->pos.vtY - 0.5);
	      iBuf[2] = 1.0 - iBuf[0];
	      iBuf[3] = 1.0 - iBuf[1];
	      rGP.dbp[0] = (vBuf[0] * iBuf[2] * iBuf[3]) +
	                   (vBuf[1] * iBuf[0] * iBuf[3]) +
			   (vBuf[2] * iBuf[2] * iBuf[1]) +
			   (vBuf[3] * iBuf[0] * iBuf[1]);
	      
	    }
	    break;
	  default:
	    break;
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  WlzGreyValueFreeWSp(gVWSp);
  if(errNum != WLZ_ERR_NONE)
  {
    if(rtnObj)
    {
      (void )WlzFreeObj(rtnObj);
      rtnObj = NULL;
    }
    else if(rtnIxv)
    {
      (void )WlzFreeIndexedValues(rtnIxv);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*!
* \return	New mesh object or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the mesh domain of the target
* 		object and (within the domain intersection) the values
* 		of the source mesh object. Given a pair of 2D conforming
* 		mesh objects. See WlzCMeshValueTransfer().
* \param	srcObj			Given source mesh object.
* \param	tgtObj			Given target mesh object.
* \param	extVal			Given external value.
* \param	interp			Interpolation method.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzCMeshValueTransferM2D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr)
{
  int		   vI,
  		   nV;
  size_t	   sz;
  WlzGreyP	   vBuf[3];
  WlzCMesh2D	   *rtnMesh,
  		   *srcMesh;
  WlzIndexedValues *rtnIxv = NULL,
  		   *srcIxv = NULL;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vBuf[0].v = vBuf[1].v = vBuf[2].v = NULL;
  srcMesh = srcObj->domain.cm2;
  rtnMesh = tgtObj->domain.cm2;
  srcIxv = srcObj->values.x;
  nV = 1;
  for(vI = 0; vI < srcIxv->rank; ++vI)
  {
    nV *= srcIxv->dim[vI];
  }
  sz = nV * sizeof(double);
  if((vBuf[0].v = AlcMalloc(sz * 3)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzValueConvertPixel(&extVal, extVal, srcIxv->vType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vBuf[1].ubp = vBuf[0].ubp + sz;
    vBuf[2].ubp = vBuf[1].ubp + sz;
    rtnIxv = WlzMakeIndexedValues(tgtObj, srcIxv->rank, srcIxv->dim,
				  srcIxv->vType, srcIxv->attach, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	rtnVal;

    rtnVal.x = rtnIxv;
    rtnObj = WlzMakeMain(tgtObj->type, tgtObj->domain, rtnVal, NULL, NULL,
                         &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcIxv->attach)
    {
      case WLZ_VALUE_ATTACH_NOD:
        {
	  int	rNI;
	  WlzCMeshEntRes *rNR,
			 *sER;

	  rNR = &(rtnMesh->res.nod);
	  sER = &(srcMesh->res.elm);
	  for(rNI = 0; rNI < rNR->maxEnt; ++rNI)
	  {
	    int		sEI;
	    WlzCMeshNod2D *rNod;

	    rNod = (WlzCMeshNod2D *)AlcVectorItemGet(rNR->vec, rNI);
	    if(rNod->idx >= 0)
	    {
	      WlzGreyP	rGP;

	      rGP.v = WlzIndexedValueGet(rtnIxv, rNI);
	      sEI = WlzCMeshElmEnclosingPos2D(srcMesh, -1,
					      rNod->pos.vtX, rNod->pos.vtY,
					      0, NULL);
	      if(sEI >= 0)
	      {
		/* Target node is in a source mesh element. */
		WlzGreyP	sGP;
		double		lambda[3];
		WlzCMeshElm2D *sElm;
		WlzCMeshNod2D *sEN[3];

		sElm = (WlzCMeshElm2D *)AlcVectorItemGet(sER->vec, sEI);
		sEN[0] =  WLZ_CMESH_ELM2D_GET_NODE_0(sElm);
		sEN[1] =  WLZ_CMESH_ELM2D_GET_NODE_1(sElm);
		sEN[2] =  WLZ_CMESH_ELM2D_GET_NODE_2(sElm);
		switch(interp)
		{
		  case WLZ_INTERPOLATION_NEAREST:
		    /* Set value to that of the closest node in the source
		     * element. */
		    {
		      int	i,
				iMin = 0;
		      double 	d,
			  	dMin;
		      WlzDVertex2 del;

		      WLZ_VTX_2_SUB(del, rNod->pos, sEN[0]->pos);
		      dMin = WLZ_VTX_2_SQRLEN(del);
		      for(i = 1; i < 3; ++i)
		      {
			WLZ_VTX_2_SUB(del, rNod->pos, sEN[i]->pos);
			d = WLZ_VTX_2_SQRLEN(del);
			if(d < dMin)
			{
			  iMin = i;
			}
		      }
		      sGP.v = WlzIndexedValueGet(srcIxv, sEN[iMin]->idx);
		      WlzValueCopyGreyToGrey(rGP, 0, srcIxv->vType,
					     sGP, 0, srcIxv->vType, nV);
		    }
		    break;
		  case WLZ_INTERPOLATION_LINEAR:
		    /* Interpolate target node value from the nodes of the
		     * source element. */
		    {
		      int		i;

		      if(WlzGeomBaryCoordsTri2D(sEN[0]->pos, sEN[1]->pos,
				    sEN[2]->pos, rNod->pos, lambda) == 0)
		      {
			lambda[0] = lambda[1] = lambda[2] = 1.0 / 3.0;
		      }
		      for(i = 0; i < 3; ++i)
		      {
			sGP.v = WlzIndexedValueGet(srcIxv, sEN[i]->idx);
			WlzValueCopyGreyToGrey(vBuf[i], 0, WLZ_GREY_DOUBLE,
					       sGP, 0, srcIxv->vType, nV);
		      }
		      for(i = 0; i < nV; ++i)
		      {
			vBuf[0].dbp[i] = (lambda[0] * vBuf[0].dbp[i]) +
				         (lambda[1] * vBuf[1].dbp[i]) +
				         (lambda[2] * vBuf[2].dbp[i]);
		      }
		      WlzValueCopyGreyToGrey(rGP, 0, srcIxv->vType,
					     vBuf[0], 0, WLZ_GREY_DOUBLE, nV);
		    }
		    break;
		  default:
		    break;
		}
	      }
	      else
	      {
		/* Target node is outside the source mesh so set the
		 * given external value. */
		WlzValueSetGrey(rGP, 0, extVal.v, srcIxv->vType, nV);
	      }
	    }
	  }
	}
	break;
      case WLZ_VALUE_ATTACH_ELM:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      default:
        errNum = WLZ_ERR_VALUES_TYPE;
	break;
    }
  }
  AlcFree(vBuf[0].v);
  if(errNum != WLZ_ERR_NONE)
  {
    if(rtnObj)
    {
      (void )WlzFreeObj(rtnObj);
      rtnObj = NULL;
    }
    else if(rtnIxv)
    {
      (void )WlzFreeIndexedValues(rtnIxv);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*!
* \return	New mesh object or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the mesh domain of the target
* 		object and (within the domain intersection) the values
* 		of the source spatial domain object. Given a pair of 3D
* 		spatial domain and conforming mesh objects. The resulting
* 		mesh values will always be WLZ_GREY_DOUBLE and attached
* 		to the mesh nodes.  See WlzCMeshValueTransfer().
* \param	srcObj			Given source spatial domain object.
* \param	tgtObj			Given target mesh object.
* \param	extVal			Given external value.
* \param	interp			Interpolation method.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzCMeshValueTransferD3D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr)
{
  WlzCMesh3D	   *rtnMesh;
  WlzIndexedValues *rtnIxv = NULL;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gVWSp->gType)
    {
      case WLZ_GREY_INT:    /* FALLTHROUGH */
      case WLZ_GREY_SHORT:  /* FALLTHROUGH */
      case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
      case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
        switch(interp)
	{
	  case WLZ_INTERPOLATION_NEAREST: /* FALLTHROUGH */
	  case WLZ_INTERPOLATION_LINEAR:
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_TYPE;
	    break;
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rtnMesh = tgtObj->domain.cm3;
    errNum = WlzValueConvertPixel(&extVal, extVal, WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rtnIxv = WlzMakeIndexedValues(tgtObj, 0, NULL,
				  WLZ_GREY_DOUBLE, WLZ_VALUE_ATTACH_NOD,
				  &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	rtnVal;

    rtnVal.x = rtnIxv;
    rtnObj = WlzMakeMain(tgtObj->type, tgtObj->domain, rtnVal, NULL, NULL,
                         &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		   rNI;
    WlzCMeshEntRes *rNR;

    rNR = &(rtnMesh->res.nod);
    for(rNI = 0; rNI < rNR->maxEnt; ++rNI)
    {
      WlzCMeshNod3D *rNod;

      rNod = (WlzCMeshNod3D *)AlcVectorItemGet(rNR->vec, rNI);
      if(rNod->idx >= 0)
      {
	WlzGreyP rGP;

	rGP.v = WlzIndexedValueGet(rtnIxv, rNI);
	switch(interp)
	{
	  case WLZ_INTERPOLATION_NEAREST:
	    WlzGreyValueGet(gVWSp,
	                    rNod->pos.vtZ, rNod->pos.vtY, rNod->pos.vtX);
	    if(gVWSp->bkdFlag)
	    {
	      rGP.dbp[0] = extVal.v.dbv;
	    }
	    else
	    {
	      switch(gVWSp->gType)
	      {
	        case WLZ_GREY_INT:
		  rGP.dbp[0] = (gVWSp->gVal[0]).inv;
		  break;
	        case WLZ_GREY_SHORT:
		  rGP.dbp[0] = (gVWSp->gVal[0]).shv;
		  break;
	        case WLZ_GREY_UBYTE:
		  rGP.dbp[0] = (gVWSp->gVal[0]).ubv;
		  break;
	        case WLZ_GREY_FLOAT:
		  rGP.dbp[0] = (gVWSp->gVal[0]).flv;
		  break;
	        case WLZ_GREY_DOUBLE:
		  rGP.dbp[0] = (gVWSp->gVal[0]).dbv;
		  break;
	        default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	    }
	    break;
	  case WLZ_INTERPOLATION_LINEAR:
	    WlzGreyValueGetCon(gVWSp,
	                       rNod->pos.vtZ, rNod->pos.vtY, rNod->pos.vtX);
	    if(gVWSp->bkdFlag)          
	    {
	      rGP.dbp[0] = extVal.v.dbv;
	    }
	    else
	    {
	      int	i;
	      double	iBuf[6],
	      		vBuf[8];

	      switch(gVWSp->gType)
	      {
	        case WLZ_GREY_INT:
		  for(i = 0; i < 8; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).inv;
		  }
		  break;
	        case WLZ_GREY_SHORT:
		  for(i = 0; i < 8; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).shv;
		  }
		  break;
	        case WLZ_GREY_UBYTE:
		  for(i = 0; i < 8; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).ubv;
		  }
		  break;
	        case WLZ_GREY_FLOAT:
		  for(i = 0; i < 8; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).flv;
		  }
		  break;
	        case WLZ_GREY_DOUBLE:
		  for(i = 0; i < 8; ++i)
		  {
		    vBuf[i] = (gVWSp->gVal[i]).dbv;
		  }
		  break;
	        default:
		  break;
	      }
	      iBuf[0] = rNod->pos.vtX - WLZ_NINT(rNod->pos.vtX - 0.5);
	      iBuf[1] = rNod->pos.vtY - WLZ_NINT(rNod->pos.vtY - 0.5);
	      iBuf[2] = rNod->pos.vtZ - WLZ_NINT(rNod->pos.vtZ - 0.5);
	      iBuf[3] = 1.0 - iBuf[0];
	      iBuf[4] = 1.0 - iBuf[1];
	      iBuf[5] = 1.0 - iBuf[2];
	      rGP.dbp[0] = (vBuf[0] * iBuf[3] * iBuf[4] * iBuf[5]) +
			   (vBuf[1] * iBuf[0] * iBuf[4] * iBuf[5]) +
			   (vBuf[2] * iBuf[3] * iBuf[1] * iBuf[5]) +
			   (vBuf[3] * iBuf[0] * iBuf[1] * iBuf[5]) +
			   (vBuf[4] * iBuf[3] * iBuf[4] * iBuf[2]) +
			   (vBuf[5] * iBuf[0] * iBuf[4] * iBuf[2]) +
			   (vBuf[6] * iBuf[3] * iBuf[1] * iBuf[2]) +
			   (vBuf[7] * iBuf[0] * iBuf[1] * iBuf[2]);
	    }
	    break;
	  default:
	    break;
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  WlzGreyValueFreeWSp(gVWSp);
  if(errNum != WLZ_ERR_NONE)
  {
    if(rtnObj)
    {
      (void )WlzFreeObj(rtnObj);
      rtnObj = NULL;
    }
    else if(rtnIxv)
    {
      (void )WlzFreeIndexedValues(rtnIxv);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*!
* \return	New mesh object or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the mesh domain of the target
* 		object and (within the domain intersection) the values
* 		of the source mesh object. Given a pair of 3D conforming
* 		mesh objects. See WlzCMeshValueTransfer().
* \param	srcObj			Given source mesh object.
* \param	tgtObj			Given target mesh object.
* \param	extVal			Given external value.
* \param	interp			Interpolation method.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzCMeshValueTransferM3D(
				  WlzObject *srcObj,
				  WlzObject *tgtObj,
				  WlzPixelV extVal,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr)
{
  int		   vI,
  		   nV;
  size_t	   sz;
  WlzGreyP	   vBuf[4];
  WlzCMesh3D	   *rtnMesh,
  		   *srcMesh;
  WlzIndexedValues *rtnIxv = NULL,
  		   *srcIxv = NULL;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vBuf[0].v = vBuf[1].v = vBuf[2].v = vBuf[3].v = NULL;
  srcMesh = srcObj->domain.cm3;
  rtnMesh = tgtObj->domain.cm3;
  srcIxv = srcObj->values.x;
  nV = 1;
  for(vI = 0; vI < srcIxv->rank; ++vI)
  {
    nV *= srcIxv->dim[vI];
  }
  sz = nV * sizeof(double);
  if((vBuf[0].v = AlcMalloc(sz * 4)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzValueConvertPixel(&extVal, extVal, srcIxv->vType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vBuf[1].ubp = vBuf[0].ubp + sz;
    vBuf[2].ubp = vBuf[1].ubp + sz;
    vBuf[3].ubp = vBuf[2].ubp + sz;
    rtnIxv = WlzMakeIndexedValues(tgtObj, srcIxv->rank, srcIxv->dim,
				  srcIxv->vType, srcIxv->attach, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	rtnVal;

    rtnVal.x = rtnIxv;
    rtnObj = WlzMakeMain(tgtObj->type, tgtObj->domain, rtnVal, NULL, NULL,
                         &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcIxv->attach)
    {
      case WLZ_VALUE_ATTACH_NOD:
        {
	  int	rNI;
	  WlzCMeshEntRes *rNR,
			 *sER;

	  rNR = &(rtnMesh->res.nod);
	  sER = &(srcMesh->res.elm);
	  for(rNI = 0; rNI < rNR->maxEnt; ++rNI)
	  {
	    int		sEI;
	    WlzCMeshNod3D *rNod;

	    rNod = (WlzCMeshNod3D *)AlcVectorItemGet(rNR->vec, rNI);
	    if(rNod->idx >= 0)
	    {
	      WlzGreyP	rGP;

	      rGP.v = WlzIndexedValueGet(rtnIxv, rNI);
	      sEI = WlzCMeshElmEnclosingPos3D(srcMesh, -1,
				  rNod->pos.vtX, rNod->pos.vtY, rNod->pos.vtZ,
				  0, NULL);
	      if(sEI >= 0)
	      {
		/* Target node is in a source mesh element. */
		WlzGreyP	sGP;
		double		lambda[4];
		WlzCMeshElm3D *sElm;
		WlzCMeshNod3D *sEN[4];

		sElm = (WlzCMeshElm3D *)AlcVectorItemGet(sER->vec, sEI);
		sEN[0] =  WLZ_CMESH_ELM3D_GET_NODE_0(sElm);
		sEN[1] =  WLZ_CMESH_ELM3D_GET_NODE_1(sElm);
		sEN[2] =  WLZ_CMESH_ELM3D_GET_NODE_2(sElm);
		sEN[3] =  WLZ_CMESH_ELM3D_GET_NODE_3(sElm);
		switch(interp)
		{
		  case WLZ_INTERPOLATION_NEAREST:
		    /* Set value to that of the closest node in the source
		     * element. */
		    {
		      int	i,
				iMin = 0;
		      double 	d,
			  	dMin;
		      WlzDVertex3 del;

		      WLZ_VTX_3_SUB(del, rNod->pos, sEN[0]->pos);
		      dMin = WLZ_VTX_3_SQRLEN(del);
		      for(i = 1; i < 4; ++i)
		      {
			WLZ_VTX_3_SUB(del, rNod->pos, sEN[i]->pos);
			d = WLZ_VTX_3_SQRLEN(del);
			if(d < dMin)
			{
			  iMin = i;
			}
		      }
		      sGP.v = WlzIndexedValueGet(srcIxv, sEN[iMin]->idx);
		      WlzValueCopyGreyToGrey(rGP, 0, srcIxv->vType,
					     sGP, 0, srcIxv->vType, nV);
		    }
		    break;
		  case WLZ_INTERPOLATION_LINEAR:
		    /* Interpolate target node value from the nodes of the
		     * source element. */
		    {
		      int		i;

		      if(WlzGeomBaryCoordsTet3D(sEN[0]->pos, sEN[1]->pos,
				                sEN[2]->pos, sEN[3]->pos,
						rNod->pos, lambda) == 0)
		      {
			lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.25;
		      }
		      for(i = 0; i < 4; ++i)
		      {
			sGP.v = WlzIndexedValueGet(srcIxv, sEN[i]->idx);
			WlzValueCopyGreyToGrey(vBuf[i], 0, WLZ_GREY_DOUBLE,
					       sGP, 0, srcIxv->vType, nV);
		      }
		      for(i = 0; i < nV; ++i)
		      {
			vBuf[0].dbp[i] = (lambda[0] * vBuf[0].dbp[i]) +
				         (lambda[1] * vBuf[1].dbp[i]) +
				         (lambda[2] * vBuf[2].dbp[i]) +
					 (lambda[3] * vBuf[3].dbp[i]);
		      }
		      WlzValueCopyGreyToGrey(rGP, 0, srcIxv->vType,
					     vBuf[0], 0, WLZ_GREY_DOUBLE, nV);
		    }
		    break;
		  default:
		    break;
		}
	      }
	      else
	      {
		/* Target node is outside the source mesh so set the
		 * given external value. */
		WlzValueSetGrey(rGP, 0, extVal.v, srcIxv->vType, nV);
	      }
	    }
	  }
	}
	break;
      case WLZ_VALUE_ATTACH_ELM:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      default:
        errNum = WLZ_ERR_VALUES_TYPE;
	break;
    }
  }
  AlcFree(vBuf[0].v);
  if(errNum != WLZ_ERR_NONE)
  {
    if(rtnObj)
    {
      (void )WlzFreeObj(rtnObj);
      rtnObj = NULL;
    }
    else if(rtnIxv)
    {
      (void )WlzFreeIndexedValues(rtnIxv);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

