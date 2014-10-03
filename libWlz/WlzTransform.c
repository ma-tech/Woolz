#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTransform_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzTransform.c
* \author       Bill Hill
* \date         April 2006
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
* \brief	Functions operating on Woolz transform unions.
* \ingroup	WlzTransform
*/

#include <Wlz.h>

/*!
* \return	New empty transform or NULL on error.
* \ingroup	WlzTransform
* \brief	Makes a new empty transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzEmptyTransform *WlzMakeEmptyTransform(WlzErrorNum *dstErr)
{
  WlzEmptyTransform *tr;
  WlzErrorNum	errNum = WLZ_ERR_MEM_ALLOC;

  if((tr = (WlzEmptyTransform *)
           AlcCalloc(1, sizeof(WlzEmptyTransform))) != NULL)
  {
    errNum = WLZ_ERR_NONE;
    tr->type = WLZ_TRANSFORM_EMPTY;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \return	New empty transform or NULL on error.
* \ingroup	WlzTransform
* \brief	Frees the given empty transform.
* \param	tr			Given empty transform.
*/
WlzErrorNum 	WlzFreeEmptyTransform(WlzEmptyTransform *tr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else if(tr->type != WLZ_TRANSFORM_EMPTY)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if(WlzUnlink(&(tr->linkcount), &errNum))
  {
    AlcFree(tr);
  }
  return(errNum);
}


/*!
* \return	Woolz error number.
* \ingroup	WlzTransform
* \brief	Free's the given Woolz transform.
* \param	tr			Given transform.
*/
WlzErrorNum	WlzFreeTransform(WlzTransform tr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr.core == NULL)
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else
  {
    switch(tr.core->type)
    {
      case WLZ_TRANSFORM_EMPTY:
	errNum = WlzFreeEmptyTransform(tr.empty);
        break;
      case WLZ_TRANSFORM_2D_AFFINE:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_REG:     /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_TRANS:   /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_NOSHEAR: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_AFFINE:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_REG:     /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_TRANS:   /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_NOSHEAR:
	errNum = WlzFreeAffineTransform(tr.affine);
        break;
      case WLZ_TRANSFORM_2D_BASISFN:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D5_BASISFN: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_BASISFN:
	errNum = WlzBasisFnFreeTransform(tr.basis);
        break;
      case WLZ_TRANSFORM_2D_MESH:	/* FALLTHROUGH */
	errNum = WlzMeshFreeTransform(tr.mesh);
        break;
      case WLZ_TRANSFORM_2D5_MESH:	/* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_MESH:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_TRANSFORM_2D_CMESH:		/* FALLTHROUGH */
      case WLZ_TRANSFORM_2D5_CMESH:		/* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_CMESH:
        errNum = WlzFreeObj(tr.obj);
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \return	New transform with type that depends on given transform
* 		types or transform with NULL core pointer on error.
* \ingroup	WlzTransform
* \brief	Computes the product of the given transforms. The resulting
* 		transform's type depends on the given transform types:
* 		\f[
		\mathbf{T_r} = \mathbf{T_0} \mathbf{T_1}
		\f]
* 		The product of any transform with an empty transform is
* 		an empty transform.
* \param	tr0			First transform (\f$\mathbf{T_0}\f$).
* \param	tr1			Second transform (\f$\mathbf{T_1}\f$).
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzTransform	WlzTransformProduct(WlzTransform tr0, WlzTransform tr1,
				    WlzErrorNum *dstErr)
{
  WlzTransform	trR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  trR.core = NULL;
  if((tr0.core == NULL) || (tr1.core == NULL))
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else
  {
    switch(tr0.core->type)                   
    {
      case WLZ_TRANSFORM_EMPTY:
        /* tr0 == WLZ_TRANSFORM_EMPTY, tr1 == WLZ_TRANSFORM_* */
	trR.empty = WlzMakeEmptyTransform(&errNum);
        break;
      case WLZ_TRANSFORM_2D_AFFINE:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_REG:     /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_TRANS:   /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_NOSHEAR: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_AFFINE:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_REG:     /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_TRANS:   /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_NOSHEAR:
	switch(tr1.core->type)                   
	{
	  case WLZ_TRANSFORM_EMPTY:
            /* tr0 == WLZ_TRANSFORM_*_AFFINE, tr1 == WLZ_TRANSFORM_EMPTY */
	    trR.empty = WlzMakeEmptyTransform(&errNum);
	    break;
	  case WLZ_TRANSFORM_2D_AFFINE:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_REG:     /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_TRANS:   /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_NOSHEAR: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_AFFINE:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_REG:     /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_TRANS:   /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_NOSHEAR:
            /* tr0 == WLZ_TRANSFORM_*_AFFINE, tr1 == WLZ_TRANSFORM_*_AFFINE */
	    trR.affine = WlzAffineTransformProduct(tr0.affine, tr1.affine,
	    				           &errNum);
	    break;
	  case WLZ_TRANSFORM_2D_BASISFN:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D5_BASISFN: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_BASISFN:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  case WLZ_TRANSFORM_2D_MESH:
            /* tr0 == WLZ_TRANSFORM_*_AFFINE, tr1 == WLZ_TRANSFORM_2D_MESH */
	    trR.mesh = WlzMeshTransformCopy(tr1.mesh, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzMeshAffineProduct(trR.mesh, tr0.affine, 0);
	    }
	    break;
	  case WLZ_TRANSFORM_2D5_MESH: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_MESH:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  case WLZ_CMESH_2D:
	  case WLZ_CMESH_2D5:
	  case WLZ_CMESH_3D:
            /* tr0 == WLZ_TRANSFORM_*_AFFINE, tr1 == WLZ_TRANSFORM_*_CMESH */
	    trR.obj = WlzCopyObject(tr1.obj, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzCMeshAffineProduct(trR.obj, tr0.affine, 1);
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_TRANSFORM_TYPE;
	    break;
	}
        break;
      case WLZ_TRANSFORM_2D_BASISFN:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D5_BASISFN: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_BASISFN:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_TRANSFORM_2D_MESH:	/* FALLTHROUGH */
	switch(tr1.core->type)                   
	{
	  case WLZ_TRANSFORM_EMPTY:
            /* tr0 == WLZ_TRANSFORM_2D_MESH, tr1 == WLZ_TRANSFORM_EMPTY */
	    trR.empty = WlzMakeEmptyTransform(&errNum);
	    break;
	  case WLZ_TRANSFORM_2D_AFFINE:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_REG:     /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_TRANS:   /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_NOSHEAR: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_AFFINE:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_REG:     /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_TRANS:   /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_NOSHEAR:
            /* tr0 == WLZ_TRANSFORM_2D_MESH, tr1 == WLZ_TRANSFORM_*_AFFINE */
	    trR.mesh = WlzMeshTransformCopy(tr0.mesh, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzMeshAffineProduct(trR.mesh, tr1.affine, 1);
	    }
	    break;
	  case WLZ_TRANSFORM_2D_BASISFN:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D5_BASISFN: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_BASISFN:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  case WLZ_TRANSFORM_2D_MESH:
            /* tr0 == WLZ_TRANSFORM_2D_MESH, tr1 == WLZ_TRANSFORM_2D_MESH */
	    trR.obj = WlzCMeshMeshMeshProduct(tr0.mesh, tr1.mesh, &errNum);
	    break;
	  case WLZ_TRANSFORM_2D5_MESH: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_MESH:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  case WLZ_CMESH_2D:
	  case WLZ_CMESH_2D5:
	  case WLZ_CMESH_3D:
            /* tr0 == WLZ_TRANSFORM_2D_MESH, tr1 == WLZ_TRANSFORM_*_CMESH */
	    trR.obj = WlzCMeshMeshProduct(tr1.obj, tr0.mesh, 1, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_TRANSFORM_TYPE;
	    break;
	}
        break;
      case WLZ_TRANSFORM_2D5_MESH:	/* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_MESH:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_2D:		/* FALLTHROUGH */
      case WLZ_CMESH_2D5:		/* FALLTHROUGH */
      case WLZ_CMESH_3D:
	switch(tr1.core->type)                   
	{
	  case WLZ_TRANSFORM_EMPTY:
            /* tr0 == WLZ_TRANSFORM_*_CMESH, tr1 == WLZ_TRANSFORM_EMPTY */
	    trR.empty = WlzMakeEmptyTransform(&errNum);
	    break;
	  case WLZ_TRANSFORM_2D_AFFINE:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_REG:     /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_TRANS:   /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D_NOSHEAR: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_AFFINE:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_REG:     /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_TRANS:   /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_NOSHEAR:
            /* tr0 == WLZ_TRANSFORM_*_CMESH, tr1 == WLZ_TRANSFORM_*_AFFINE */
	    trR.obj = WlzCopyObject(tr0.obj, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzCMeshAffineProduct(trR.obj, tr1.affine, 0);
	    }
	    break;
	  case WLZ_TRANSFORM_2D_BASISFN:  /* FALLTHROUGH */
	  case WLZ_TRANSFORM_2D5_BASISFN: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_BASISFN:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  case WLZ_TRANSFORM_2D_MESH:
            /* tr0 == WLZ_TRANSFORM_*_CMESH, tr1 == WLZ_TRANSFORM_2D_MESH */
	    trR.obj = WlzCMeshMeshProduct(tr0.obj, tr1.mesh, 0, &errNum);
	    break;
	  case WLZ_TRANSFORM_2D5_MESH: /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_MESH:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  case WLZ_CMESH_2D:
	  case WLZ_CMESH_2D5:
	  case WLZ_CMESH_3D:
            /* tr0 == WLZ_TRANSFORM_*_CMESH, tr1 == WLZ_TRANSFORM_*_CMESH */
	    trR.obj = WlzCMeshProduct(tr0.obj, tr1.obj, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_TRANSFORM_TYPE;
	    break;
	}
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trR);
}
