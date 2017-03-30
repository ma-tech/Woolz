#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshDispField_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzCMeshDispField.c
* \author       Bill Hill
* \date         September 2016
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
* \brief	Functions transferring displacements between meshes
* 		and displacement fields.
* \ingroup	WlzMesh
*/

#include <limits.h>
#include <float.h>
#include <math.h>
#include <Wlz.h>

/*!
* \return	New constrained mesh object with displacements set
* 		from the displacement field object or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new 2 or 3D constrained mesh object using the
* 		domain of the given mesh object (with possible refinement)
* 		and creating new values for the domain which are set from
* 		the displacement field.
* 		Mesh refinement will be done if the fractional error
* 		length at element centroids is greater than the given
* 		maximum and the resulting element will not have an edge
* 		length less than the given minimum and the given element
* 		does not have a maximum to minimum edge length ratio greater
* 		than two.
* \param	mObj			Given constrained mesh object.
* \param	dObj			Given displacement field object
* 					which should be a compound array
* 					object with the appropriate number
* 					of components for the mesh dimension.
* \param	bgd			Background displacement for mesh
* 					nodes which are outside the given
* 					displacement field object.
* \param	itp			Interpolation method (only
* 					WLZ_INTERPOLATION_NEAREST and
* 					WLZ_INTERPOLATION_LINEAR are valid).
* \param	abs			Non-zero if the displacements are
* 					absolute positions rather than relative
* 					displacements.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzCMeshSetDispFromField(
				  WlzObject *mObj,
				  WlzObject *dObj,
				  WlzDVertex3 bgd,
				  WlzInterpolationType itp,
				  int abs,
				  WlzErrorNum *dstErr)
{

  int		idN,
  		dim = 0;
  WlzObject	*rObj = NULL;
  WlzIndexedValues *ixv = NULL;
  WlzCompoundArray *dCpd = NULL;
  WlzGreyValueWSpace *gVWSp[3] = {NULL};
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((mObj == NULL) || (dObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((itp != WLZ_INTERPOLATION_NEAREST) && 
          (itp != WLZ_INTERPOLATION_LINEAR))
  {
    errNum = WLZ_ERR_INTERPOLATION_TYPE;
  }
  else if(((mObj->type != WLZ_CMESH_2D) && (mObj->type != WLZ_CMESH_3D)) ||
          (dCpd = (WlzCompoundArray *)dObj)->type != WLZ_COMPOUND_ARR_1)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else
  {
    WlzObjectType dType;

    if(mObj->type == WLZ_CMESH_2D)
    {
      dim = 2;
      dType = WLZ_2D_DOMAINOBJ;
    }
    else /* mObj->type == WLZ_CMESH_3D */
    {
      dim = 3;
      dType = WLZ_3D_DOMAINOBJ;
    }
    if(dCpd->n < dim)
    {
      errNum = WLZ_ERR_OBJECT_DATA;
    }
    for(idN = 0; (errNum == WLZ_ERR_NONE) && (idN < dim); ++idN)
    {
      if(dCpd->o[idN] == NULL)
      {
	errNum = WLZ_ERR_OBJECT_NULL;
      }
      else if(dCpd->o[idN]->type != dType)
      {
	errNum = WLZ_ERR_OBJECT_TYPE;
      }
      else if(dCpd->o[idN]->domain.core == NULL)
      {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if(dCpd->o[idN]->values.core == NULL)
      {
	errNum = WLZ_ERR_VALUES_NULL;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	gVWSp[idN] = WlzGreyValueMakeWSp(dCpd->o[idN], &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(idN == 0)
	{
	  switch(gVWSp[0]->gType)
	  {
	    case WLZ_GREY_INT:    /* FALLTHROUGH */
	    case WLZ_GREY_SHORT:  /* FALLTHROUGH */
	    case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
	    case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
	    case WLZ_GREY_DOUBLE: /* FALLTHROUGH */
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
        else
	{
	  if((gVWSp[idN]->gType != gVWSp[idN - 1]->gType))
	  {
	    errNum = WLZ_ERR_GREY_TYPE;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ixv = WlzMakeIndexedValues(mObj, 1, &dim, WLZ_GREY_DOUBLE,
                               WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues val;

    val.x = ixv;
    rObj = WlzMakeMain(mObj->type, mObj->domain, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshEntRes *mnr;

    if(dim == 2)
    {
      WlzCMesh2D *mesh;

      mesh = mObj->domain.cm2;
      mnr = &(mesh->res.nod);
      for(idN = 0; idN < mnr->maxEnt; ++idN)
      {
        WlzCMeshNod2D	*nd;

        nd = (WlzCMeshNod2D *)AlcVectorItemGet(mnr->vec, idN);
	if(nd->idx >= 0)
	{
	  double 	*v;

	  v = (double *)WlzIndexedValueGet(ixv, idN);
	  if(itp == WLZ_INTERPOLATION_NEAREST)
	  {
	    int		idC;

	    for(idC = 0; idC < 2; ++idC)
	    {
	      WlzGreyValueGet(gVWSp[idC],
	                      0.0, nd->pos.vtY, nd->pos.vtX);
	    }
	    switch(gVWSp[0]->gType)
	    {
	      case WLZ_GREY_INT:
		for(idC = 0; idC < 2; ++idC)
		{
	          v[idC] = gVWSp[idC]->gVal[0].inv;
	        }
		break;
	      case WLZ_GREY_SHORT:
		for(idC = 0; idC < 2; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].shv;
		}
		break;
	      case WLZ_GREY_UBYTE:
		for(idC = 0; idC < 2; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].ubv;
		}
		break;
	      case WLZ_GREY_FLOAT:
		for(idC = 0; idC < 2; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].flv;
		}
		break;
	      case WLZ_GREY_DOUBLE:
		for(idC = 0; idC < 2; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].dbv;
		}
	        break;
	      default:
	        break;
	    }
	  }
	  else /* itp == WLZ_INTERPOLATION_LINEAR */
	  {
	    int		idC,
			idD;
	    WlzDVertex2 p0,
			p1;
	    double	q[4];

	    for(idC = 0; idC < 3; ++idC)
	    {
	      v[idC] = 0.0;
	      WlzGreyValueGetCon(gVWSp[idC],
	                         0.0, nd->pos.vtY, nd->pos.vtX);
	    }
	    p0.vtX = nd->pos.vtX - WLZ_NINT(nd->pos.vtX - 0.5);
	    p0.vtY = nd->pos.vtY - WLZ_NINT(nd->pos.vtY - 0.5);
	    p1.vtX = 1.0 - p0.vtX;
	    p1.vtY = 1.0 - p0.vtY;
	    q[0] = p1.vtX * p1.vtY;
	    q[1] = p1.vtX * p0.vtY;
	    q[2] = p1.vtX * p0.vtY;
	    q[3] = p0.vtX * p0.vtY;
	    switch(gVWSp[0]->gType)
	    {
	      case WLZ_GREY_INT:
		for(idC = 0; idC < 2; ++idC)
		{
		  for(idD = 0; idD < 4; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].inv;
		  }
		}
		break;
	      case WLZ_GREY_SHORT:
		for(idC = 0; idC < 2; ++idC)
		{
		  for(idD = 0; idD < 4; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].shv;
		  }
		}
		break;
	      case WLZ_GREY_UBYTE:
		for(idC = 0; idC < 2; ++idC)
		{
		  for(idD = 0; idD < 4; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].ubv;
		  }
		}
		break;
	      case WLZ_GREY_FLOAT:
		for(idC = 0; idC < 2; ++idC)
		{
		  for(idD = 0; idD < 4; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].flv;
		  }
		}
		break;
	      case WLZ_GREY_DOUBLE:
		for(idC = 0; idC < 2; ++idC)
		{
		  for(idD = 0; idD < 4; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].dbv;
		  }
		}
		break;
	      default:
		break;
	    }
	  }
	  if(!abs)
	  {
	    v[0] -= nd->pos.vtX;
	    v[1] -= nd->pos.vtY;
	  }
        }
      }
    }
    else /* dim == 3 */
    {
      WlzCMesh3D *mesh;

      mesh = mObj->domain.cm3;
      mnr = &(mesh->res.nod);
      for(idN = 0; idN < mnr->maxEnt; ++idN)
      {
        WlzCMeshNod3D	*nd;

        nd = (WlzCMeshNod3D *)AlcVectorItemGet(mnr->vec, idN);
	if(nd->idx >= 0)
	{
	  double 	*v;

	  v = (double *)WlzIndexedValueGet(ixv, idN);
	  if(itp == WLZ_INTERPOLATION_NEAREST)
	  {
	    int		idC;

	    for(idC = 0; idC < 3; ++idC)
	    {
	      WlzGreyValueGet(gVWSp[idC],
	                      nd->pos.vtZ, nd->pos.vtY, nd->pos.vtX);
	    }
	    switch(gVWSp[0]->gType)
	    {
	      case WLZ_GREY_INT:
		for(idC = 0; idC < 3; ++idC)
		{
	          v[idC] = gVWSp[idC]->gVal[0].inv;
	        }
		break;
	      case WLZ_GREY_SHORT:
		for(idC = 0; idC < 3; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].shv;
		}
		break;
	      case WLZ_GREY_UBYTE:
		for(idC = 0; idC < 3; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].ubv;
		}
		break;
	      case WLZ_GREY_FLOAT:
		for(idC = 0; idC < 3; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].flv;
		}
		break;
	      case WLZ_GREY_DOUBLE:
		for(idC = 0; idC < 3; ++idC)
		{
		  v[idC] = gVWSp[idC]->gVal[0].dbv;
		}
	        break;
	      default:
	        break;
	    }
	  }
	  else /* itp == WLZ_INTERPOLATION_LINEAR */
	  {
	    int		idC,
			idD;
	    WlzDVertex3 p0,
			p1;
	    double	q[8];

	    for(idC = 0; idC < 3; ++idC)
	    {
	      v[idC] = 0.0;
	      WlzGreyValueGetCon(gVWSp[idC],
	                         nd->pos.vtZ, nd->pos.vtY, nd->pos.vtX);
	    }
	    p0.vtX = nd->pos.vtX - WLZ_NINT(nd->pos.vtX - 0.5);
	    p0.vtY = nd->pos.vtY - WLZ_NINT(nd->pos.vtY - 0.5);
	    p0.vtZ = nd->pos.vtZ - WLZ_NINT(nd->pos.vtZ - 0.5);
	    p1.vtX = 1.0 - p0.vtX;
	    p1.vtY = 1.0 - p0.vtY;
	    p1.vtZ = 1.0 - p0.vtZ;
	    q[0] = p1.vtX * p1.vtY * p1.vtZ;
	    q[1] = p0.vtX * p1.vtY * p1.vtZ;
	    q[2] = p1.vtX * p0.vtY * p1.vtZ;
	    q[3] = p0.vtX * p0.vtY * p1.vtZ;
	    q[4] = p1.vtX * p1.vtY * p0.vtZ;
	    q[5] = p0.vtX * p1.vtY * p0.vtZ;
	    q[6] = p1.vtX * p0.vtY * p0.vtZ;
	    q[7] = p0.vtX * p0.vtY * p0.vtZ;
	    switch(gVWSp[0]->gType)
	    {
	      case WLZ_GREY_INT:
		for(idC = 0; idC < 3; ++idC)
		{
		  for(idD = 0; idD < 8; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].inv;
		  }
		}
		break;
	      case WLZ_GREY_SHORT:
		for(idC = 0; idC < 3; ++idC)
		{
		  for(idD = 0; idD < 8; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].shv;
		  }
		}
		break;
	      case WLZ_GREY_UBYTE:
		for(idC = 0; idC < 3; ++idC)
		{
		  for(idD = 0; idD < 8; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].ubv;
		  }
		}
		break;
	      case WLZ_GREY_FLOAT:
		for(idC = 0; idC < 3; ++idC)
		{
		  for(idD = 0; idD < 8; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].flv;
		  }
		}
		break;
	      case WLZ_GREY_DOUBLE:
		for(idC = 0; idC < 3; ++idC)
		{
		  for(idD = 0; idD < 8; ++idD)
		  {
		    v[idC] += q[idD] * gVWSp[idC]->gVal[idD].dbv;
		  }
		}
		break;
	      default:
		break;
	    }
	  }
	  if(!abs)
	  {
	    v[0] -= nd->pos.vtX;
	    v[1] -= nd->pos.vtY;
	    v[2] -= nd->pos.vtZ;
	  }
        }
      }
    }
  }
  for(idN = 0; idN < dim; ++idN)
  {
    WlzGreyValueFreeWSp(gVWSp[idN]);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}


/*!
* \return	New displacement field object (a compound array object)
* 		or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new 2 or 3D compound array object with the
* 		component displacements of the given mesh interpolated
* 		over the mesh domain.
* \param	mObj			Given constrained mesh object
* 					with displacement values.
* \param	bgd			Background displacement for field
* 					values outside the given mesh.
* \param	itp			Interpolation method (only
* 					WLZ_INTERPOLATION_NEAREST and
* 					WLZ_INTERPOLATION_LINEAR are valid).
* \param 	invert			Invert the transform if non-zero.
* \param	abs			Non-zero if the field displacements
* 					absolute positions rather than the
* 					relative displacements in the mesh
* 					values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzCMeshDispToField(
				  WlzObject *mObj,
				  WlzDVertex3 bgd,
				  WlzInterpolationType itp,
				  int invert,
				  int abs,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((itp != WLZ_INTERPOLATION_NEAREST) &&
          (itp != WLZ_INTERPOLATION_LINEAR))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((mObj->values.core->type != (WlzObjectType )WLZ_INDEXED_VALUES) ||
          ((mObj->values.x->attach != WLZ_VALUE_ATTACH_NOD) &&
	   (mObj->values.x->attach != WLZ_VALUE_ATTACH_ELM)))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    switch(mObj->type)
    {
      case WLZ_CMESH_2D: /* FALLTHROUGH */
      case WLZ_CMESH_3D:
	{
	  int		dim;
	  WlzObjectType   oType;
	  WlzObject	*dObj = NULL,
			  *nmObj = NULL;
	  WlzCompoundArray *rCpd = NULL;
	  WlzErrorNum	errNum = WLZ_ERR_NONE;

	  if(mObj->type == WLZ_CMESH_2D)
	  {
	    dim = 2;
	    oType = WLZ_2D_DOMAINOBJ;
	  }
	  else
	  {
	    dim = 3;
	    oType = WLZ_3D_DOMAINOBJ;
	  }
	  dObj = WlzAssignObject(
		 WlzCMeshToDomObj(mObj, 0, 1.0, &errNum), NULL);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    rCpd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, dim, NULL,
					oType, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzValues	nVal;

	    if(!abs)
	    {
	      nVal = mObj->values;
	    }
	    else
	    {
	      int		maxEnt = 0;
	      size_t	vCnt = 0;
	      AlcVector	*eVec = NULL,
			  *vVec0 = NULL,
			  *vVec1 = NULL;

	      (void )WlzIndexedValueSize(mObj->values.x, &vCnt, &errNum);
	      if(errNum == WLZ_ERR_NONE) 
	      {
		nVal.x = WlzCopyIndexedValues(mObj->values.x, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE) 
	      {
		int	ide;

		vVec0 = mObj->values.x->values;
		vVec1 = nVal.x->values;
		if(dim == 2)
		{
		  switch(mObj->values.x->attach)
		  {
		    case WLZ_VALUE_ATTACH_NOD:
		      eVec = mObj->domain.cm2->res.nod.vec;
		      maxEnt = mObj->domain.cm2->res.nod.maxEnt;
		      break;
		    case WLZ_VALUE_ATTACH_ELM:
		      eVec = mObj->domain.cm2->res.elm.vec;
		      maxEnt = mObj->domain.cm2->res.elm.maxEnt;
		      break;
		    default:
		      break;
		  }
		}
		else /* dim == 3 */
		{
		  switch(mObj->values.x->attach)
		  {
		    case WLZ_VALUE_ATTACH_NOD:
		      eVec = mObj->domain.cm3->res.nod.vec;
		      maxEnt = mObj->domain.cm3->res.nod.maxEnt;
		      break;
		    case WLZ_VALUE_ATTACH_ELM:
		      eVec = mObj->domain.cm3->res.elm.vec;
		      maxEnt = mObj->domain.cm3->res.elm.maxEnt;
		      break;
		    default:
		      break;
		  }
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(ide = 0; ide < maxEnt; ++ide)
		{
		  WlzCMeshEntP entP;

		  if(errNum == WLZ_ERR_NONE)
		  {
		    entP.v = AlcVectorItemGet(eVec, ide);
		    if(entP.core->idx >= 0)
		    {
		      int		idv;
		      WlzGreyP 	p0,
				  p1;

		      p0.v = AlcVectorItemGet(vVec0, ide);
		      p1.v = AlcVectorItemGet(vVec1, ide);
		      switch(mObj->values.x->vType)
		      {
			case WLZ_GREY_INT:
			  for(idv = 0; idv < vCnt; ++idv)
			  {
			    p1.inp[idv] -= p0.inp[idv];
			  }
			  break;
			case WLZ_GREY_SHORT:
			  for(idv = 0; idv < vCnt; ++idv)
			  {
			    p1.shp[idv] -= p0.shp[idv];
			  }
			  break;
			case WLZ_GREY_UBYTE:
			  for(idv = 0; idv < vCnt; ++idv)
			  {
			    p1.ubp[idv] -= p0.ubp[idv];
			  }
			  break;
			case WLZ_GREY_FLOAT:
			  for(idv = 0; idv < vCnt; ++idv)
			  {
			    p1.flp[idv] -= p0.flp[idv];
			  }
			  break;
			case WLZ_GREY_DOUBLE:
			  for(idv = 0; idv < vCnt; ++idv)
			  {
			    p1.dbp[idv] -= p0.dbp[idv];
			  }
			  break;
			default:
#ifdef _OPENMP
#pragma omp critical (WlzCMeshDispToField)
			  {
			    errNum = WLZ_ERR_GREY_TYPE;
			  }
#else
                          errNum = WLZ_ERR_GREY_TYPE;
#endif
			  break;
		      }
		    }
		  }
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      nmObj = WlzMakeMain(mObj->type, mObj->domain, nVal, NULL, NULL,
				  &errNum);
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int 	idC;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for(idC = 0; idC < dim; ++idC)
	    {
	      if(errNum == WLZ_ERR_NONE)
	      {
		WlzErrorNum errNum2 = WLZ_ERR_NONE;

		rCpd->o[idC] = WlzAssignObject(
			       WlzCMeshToDomObjValues(dObj, nmObj, itp, idC,
						      &errNum2), NULL);
		if(errNum2 == WLZ_ERR_NONE)
		{
		  WlzPixelV bgdV;

		  bgdV.type = WLZ_GREY_DOUBLE;
		  switch(idC)
		  {
		    case 0:
		      bgdV.v.dbv = bgd.vtX;
		      break;
		    case 1:
		      bgdV.v.dbv = bgd.vtY;
		      break;
		    case 2:
		      bgdV.v.dbv = bgd.vtZ;
		      break;
		    default:
		      break;
		  }
		  errNum2 = WlzSetBackground(rCpd->o[idC], bgdV);
		}
		if(errNum2 != WLZ_ERR_NONE)
		{
#ifdef _OPENMP
#pragma omp critical (WlzCMeshDispToField)
		  {
		    if(errNum == WLZ_ERR_NONE)
		    {
		      errNum = errNum2;
		    }
		  }
#else
		  errNum = errNum2;
#endif
		}
	      }
	    }
	  }
	  (void )WlzFreeObj(dObj);
	  (void )WlzFreeObj(nmObj);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    (void )WlzFreeObj((WlzObject *)rCpd);
	    rCpd = NULL;
	  }
	  if(dstErr)
	  {
	    *dstErr = errNum;
	  }
	  rObj = (WlzObject *)rCpd;
	}
	break;
      case WLZ_CMESH_2D5:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}
