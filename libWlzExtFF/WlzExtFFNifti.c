#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFAnl_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFAnl.c
* \author       Bill Hill
* \date         February 2005
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
* \brief	Functions for reading and writting Woolz objects to
* 		the ANALYZE 7.5 format file.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>
#include <nifti1_io.h>

static void			WlzEffNiftiCopyToObj1D(
				  WlzGreyP wGP,
				  WlzGreyP nGP,
				  int nDT,
				  int nVPP,
				  int nBPP,
				  int nV,
				  int idx);
static char 			*WlzEffNiftiInfoToStr(
				  nifti_image *nim,
				  WlzErrorNum *dstErr);
static int 			WlzEffNiftiFromWlzGType(
				  WlzGreyType wGType,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzEffNiftiToObj3D(
				  nifti_image *nim,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzEffNiftiToObj2D(
				  nifti_image *nim,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzEffNiftiCopyToObj3D(
				  WlzGreyP nGP,
				  int nDT,
				  int nBPP,
				  WlzGreyType wGType,
				  int nVPP,
				  WlzPixelV bgdV,
				  WlzIVertex3 sz,
				  int idx,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzEffNiftiCopyToObj2D(
				  WlzGreyP nGP,
				  int nDT,
				  int nBPP,
				  WlzGreyP wGP,
				  WlzGreyType wGType,
				  int nVPP,
				  WlzPixelV bgdV,
				  WlzIVertex2 sz,
				  int idx,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzEffNiftiToWlzType(
				  nifti_image *nim,
				  int *dstNVPP,
				  WlzGreyType *dstWGType);
/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
*		NIfTI format. In some cases (for example when the NIfTI
*		object has complex values) a compound object may be
*		returned. A text property is created for the returned
*		object, containing various fields of the NIfTI header,
*		some of which may not be used in the conversion. The
*		qform transform/quaternion, mapping the image to it's
*		source (scanner) coordinates is copied text property
*		but no transform is either created of applied for it.
*		The sform transform/quaternion is always applied if the
*		parameter flag is set, otherwise it ionly is applied if
*		it is a pure translation. If the sform transform is not
*		applied then it will (if not an identity transform) be
*		used to create a trans-obj. If affine transforms are
*		applied then WLZ_INTERPOLATION_NEAREST is used. NIfTI
*		allows the linear transformation of grey values. This
*		will only be applied if the parameter flag is set. The
*		grey scaling will result in floating point image grey
*		values.
* \param	gvnFileName		Given file name.
* \param	spatialTr		Apply the NIfTI sform spatial
* 					transform (mapping the image to some
* 					global coordinate system) to the image.
* 					This is not required to recover a
* 					simple offset.
* \param	greySc			Apply the NIfTI grey scaling.
* \param	dstErr			Destination error code, may be NULL.
*/
WlzObject	*WlzEffReadObjNifti(const char *gvnFileName,
				    int spatialTr,
				    int greySc,
				    WlzErrorNum *dstErr)
{
  nifti_image	*nim = NULL;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	trEps = 0.000001;

  /* Read NIfTI object into NIfTI image format. */
  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    nim = nifti_image_read(gvnFileName, 1);
    if(nim == NULL)
    {
      errNum = WLZ_ERR_READ_EOF;
    }
  }
  /* Create basic Woolz domain object from the NIfTI image. */
  if(errNum == WLZ_ERR_NONE)
  {
    switch(nim->ndim)
    {
      case 2:
        obj = WlzEffNiftiToObj2D(nim, &errNum);
	break;
      case 3:
        obj = WlzEffNiftiToObj3D(nim, &errNum);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Apply grey value scale and offset if required. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((greySc != 0) &&
       (fabs(nim->scl_slope) > FLT_MIN) &&
       ((fabs(nim->scl_slope - 1.0) > FLT_MIN) ||
        (fabs(nim->scl_inter) > FLT_MIN)))
    {
      WlzPixelV	a,
      		m;
      WlzObject	*nObj;

      a.type = m.type = WLZ_GREY_FLOAT;
      a.v.flv = nim->scl_inter;
      m.v.flv = nim->scl_slope;
      nObj = WlzScalarMulAdd(obj, m, a, WLZ_GREY_FLOAT, &errNum);
      WlzFreeObj(obj);
      obj = nObj;
    }
  }
  /* Create affine transform for the s forms if needed.
   * q is from object to scanner space which is put into the NIfTI
   * info text property, s is from object to some standard space
   * and is encoded in a WlzAffineTransform.. */
  if((errNum == WLZ_ERR_NONE) && (nim->sform_code > 0))
  {
    int		idI,
		idJ;
    WlzAffineTransform *sTr = NULL;
    WlzTransformType sTrType = WLZ_TRANSFORM_3D_AFFINE;

    sTr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      for(idJ = 0; idJ < 4; ++idJ)
      {
	for(idI = 0; idI < 4; ++idI)
	{
	  sTr->mat[idJ][idI] = nim->sto_ijk.m[idJ][idI];
	}
      }
      if((fabs(nim->sto_ijk.m[0][0] - 1.0) < trEps) &&
	 (fabs(nim->sto_ijk.m[1][1] - 1.0) < trEps) &&
	 (fabs(nim->sto_ijk.m[2][2] - 1.0) < trEps) &&
	 (fabs(nim->sto_ijk.m[3][3] - 1.0) < trEps) &&
	 (fabs(nim->sto_ijk.m[0][1]) < trEps) &&
	 (fabs(nim->sto_ijk.m[0][2]) < trEps) &&
	 (fabs(nim->sto_ijk.m[1][0]) < trEps) &&
	 (fabs(nim->sto_ijk.m[1][2]) < trEps) &&
	 (fabs(nim->sto_ijk.m[2][0]) < trEps) &&
	 (fabs(nim->sto_ijk.m[2][1]) < trEps) &&
	 (fabs(nim->sto_ijk.m[3][0]) < trEps) &&
	 (fabs(nim->sto_ijk.m[3][1]) < trEps) &&
	 (fabs(nim->sto_ijk.m[3][2]) < trEps))
      {
	if((nim->ndim == 2) && (fabs(nim->sto_ijk.m[2][3]) < trEps))
	{
	  if((fabs(nim->sto_ijk.m[0][3]) < trEps) &&
	     (fabs(nim->sto_ijk.m[1][3]) < trEps))
	  {
	    sTrType = WLZ_TRANSFORM_EMPTY;
	  }
	  else
	  {
	    sTrType = WLZ_TRANSFORM_2D_TRANS;
	  }
	}
	else
	{
	  if((fabs(nim->sto_ijk.m[0][3]) < trEps) &&
	     (fabs(nim->sto_ijk.m[1][3]) < trEps) &&
	     (fabs(nim->sto_ijk.m[2][3]) < trEps))
	  {
	    sTrType = WLZ_TRANSFORM_EMPTY;
	  }
	  else if(nim->ndim == 3)
	  {
	    sTrType = WLZ_TRANSFORM_3D_TRANS;
	  }
	}
      }
      if(((nim->ndim == 2) && (sTrType == WLZ_TRANSFORM_2D_TRANS)) ||
         ((nim->ndim == 3) && (sTrType == WLZ_TRANSFORM_3D_TRANS)))
      {
        /* Shift the 2D object(s). */
	int	  idN;
	WlzIVertex3 sft;
	WlzObject *nObj;

	if(nim->ndim == 2)
	{
          sft.vtX = sTr->mat[0][2];
	  sft.vtY = sTr->mat[1][2];
	  sft.vtZ = 0;
	}
	else /* nim->ndim == 3 */
	{
          sft.vtX = sTr->mat[0][3];
	  sft.vtY = sTr->mat[1][3];
	  sft.vtZ = sTr->mat[2][3];
	}
        switch(obj->type)
	{
	  case WLZ_COMPOUND_ARR_1:
	    for(idN = 0; idN < ((WlzCompoundArray *)obj)->n; ++idN)
	    {
	      nObj = WlzShiftObject(((WlzCompoundArray *)obj)->o[idN],
	                            sft.vtX, sft.vtY, sft.vtZ,
				    &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		(void )WlzFreeObj(((WlzCompoundArray *)obj)->o[idN]);
		((WlzCompoundArray *) obj)->o[idN] =
				WlzAssignObject(nObj, NULL);
	      }
	      else
	      {
	        break;
	      }
	    }
	  case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
	  case WLZ_3D_DOMAINOBJ:
	    nObj = WlzShiftObject(obj, sft.vtX, sft.vtY, sft.vtZ,
	                          &errNum);
	    
	    if(errNum == WLZ_ERR_NONE)
	    {
	      (void )WlzFreeObj(obj);
	      obj = nObj;
	    }
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
      else if(sTrType == WLZ_TRANSFORM_3D_AFFINE)
      {
        if(spatialTr == 0)
	{
	  /* Create a trans-obj. */
	  WlzProperty p;
	  WlzPropertyList *pLst = NULL;

	  p.name = WlzMakeNameProperty("NIfTI sto_ijk", &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if((pLst = WlzMakePropertyList(NULL)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(AlcDLPListEntryAppend(pLst->list, NULL, (void *)(p.core),
				     WlzFreePropertyListEntry) != ALC_ER_NONE)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzDomain dom;
	    WlzValues val;
	    WlzObject	*nObj;

	    dom.t = sTr;
	    val.obj = obj;
	    nObj = WlzMakeMain(WLZ_TRANS_OBJ, dom, val, pLst, NULL, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      obj = nObj;
	      sTr = NULL;
	    }
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    if(pLst != NULL)
	    {
	      WlzFreePropertyListEntry(pLst);
	    }
	    else if(p.core != NULL)
	    {
	      (void )WlzFreeProperty(p);
	    }
	  }
	}
	else
	{
          /* Affine transform the 2 or 3D object(s). */
	  int	  idN;
	  WlzObject *nObj;

	  switch(obj->type)
	  {
	    case WLZ_COMPOUND_ARR_1:
	      for(idN = 0; idN < ((WlzCompoundArray *)obj)->n; ++idN)
	      {
		nObj = WlzAffineTransformObj(((WlzCompoundArray *)obj)->o[idN],
				      sTr, WLZ_INTERPOLATION_NEAREST,
				      &errNum);
		if(errNum == WLZ_ERR_NONE)
		{
		  (void )WlzFreeObj(((WlzCompoundArray *)obj)->o[idN]);
		  ((WlzCompoundArray *) obj)->o[idN] =
				  WlzAssignObject(nObj, NULL);
		}
		else
		{
		  break;
		}
	      }
	    case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
	    case WLZ_3D_DOMAINOBJ:
	      nObj = WlzAffineTransformObj(((WlzCompoundArray *)obj)->o[idN],
				    sTr, WLZ_INTERPOLATION_NEAREST,
				    &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		(void )WlzFreeObj(obj);
		obj = nObj;
	      }
	    default:
	      errNum = WLZ_ERR_OBJECT_TYPE;
	      break;
	  }
	}
      }
    }
    (void )WlzFreeAffineTransform(sTr);
  }
  /* Create properties with intent code and descriptive text. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzProperty	p;
    char 	*s = NULL;
    WlzPropertyList *pLst;

    p.core = NULL;
    nim->descrip[79] = '\0';
    nim->intent_name[15] = '\0';
    s = WlzEffNiftiInfoToStr(nim, &errNum);
    if((errNum == WLZ_ERR_NONE) &&
       ((pLst = WlzMakePropertyList(NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      p.text = WlzMakeTextProperty("NIfTI Info", s, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)             
    {
      if(AlcDLPListEntryAppend(pLst->list, NULL, (void *)(p.core),
      			       WlzFreePropertyListEntry) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      obj->plist = WlzAssignPropertyList(pLst, NULL);
    }
    else
    {
      if(pLst != NULL)
      {
        (void )WlzFreePropertyList(pLst);
      }
    }
    AlcFree(s);
  }
  if(nim != NULL)
  {
    nifti_image_free(nim);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file(s)
*		using the NIfTI file format.
* \param	gvnFileName		Given file name with .hdr, .img or no
* 					extension.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjNifti(const char *gvnFileName, WlzObject *obj)
{
  int		nDType;
  nifti_image	*nim = NULL;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gType = WlzGreyTypeFromObj(obj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nDType = WlzEffNiftiFromWlzGType(gType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	{
	  WlzIBox2 bBox;
	  void	  **datAry = NULL; 

	  bBox.xMin = obj->domain.i->kol1;
	  bBox.xMax = obj->domain.i->lastkl;
	  bBox.yMin = obj->domain.i->line1;
	  bBox.yMax = obj->domain.i->lastln;
	  if((nim = (nifti_image *)AlcCalloc(1, sizeof(nifti_image))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
            nim->dim[0] = nim->ndim = 2;
	    nim->dim[1] = nim->nx = bBox.xMax - bBox.xMin + 1;
	    nim->dim[2] = nim->ny = bBox.yMax - bBox.yMin + 1;
	    nim->dim[3] = nim->nz = 1;
	    nim->dim[4] = nim->nt = 1;
	    nim->dim[5] = nim->nu = 1;
	    nim->dim[6] = nim->nv = 1;
	    nim->dim[7] = nim->nw = 1;
	    nim->nvox = nim->nx * nim->ny;
	    nim->datatype = nDType;
	    nim->pixdim[1] = nim->dx = 1.0;
	    nim->pixdim[2] = nim->dy = 1.0;
	    nim->pixdim[3] = nim->dz = 1.0;
	    nim->pixdim[4] = nim->dt = 1.0;
	    nim->pixdim[5] = nim->du = 1.0;
	    nim->pixdim[6] = nim->dv = 1.0;
	    nim->pixdim[7] = nim->dw = 1.0;
	    nim->scl_slope = 1.0;
	    nim->scl_inter = 0.0;
	    nim->sform_code = 1;
	    nim->qto_xyz.m[0][0] = 1.0;
	    nim->qto_xyz.m[1][1] = 1.0;
	    nim->qto_xyz.m[2][2] = 1.0;
	    nim->qto_xyz.m[3][3] = 1.0;
	    nim->qto_ijk.m[0][0] = 1.0;
	    nim->qto_ijk.m[1][1] = 1.0;
	    nim->qto_ijk.m[2][2] = 1.0;
	    nim->qto_ijk.m[3][3] = 1.0;
	    nim->sto_xyz.m[0][0] = 1.0; nim->sto_xyz.m[0][3] = -(bBox.xMin);
	    nim->sto_xyz.m[1][1] = 1.0; nim->sto_xyz.m[1][3] = -(bBox.yMin);
	    nim->sto_xyz.m[2][2] = 1.0;
	    nim->sto_xyz.m[3][3] = 1.0;
	    nim->sto_ijk.m[0][0] = 1.0; nim->sto_ijk.m[0][3] = bBox.xMin;
	    nim->sto_ijk.m[1][1] = 1.0; nim->sto_ijk.m[1][3] = bBox.yMin;
	    nim->sto_ijk.m[2][2] = 1.0;
	    nim->sto_ijk.m[3][3] = 1.0;
	    nifti_mat44_to_quatern(nim->qto_xyz,
	                           &(nim->quatern_b),
	                           &(nim->quatern_c),
	                           &(nim->quatern_d),
	                           &(nim->qoffset_x),
	                           &(nim->qoffset_y),
	                           &(nim->qoffset_z),
				   NULL, NULL, NULL, &(nim->qfac));
	    nifti_datatype_sizes(nim->datatype,
	                         &(nim->nbyper), &(nim->swapsize) ) ;
            nim->byteorder = nifti_short_order();
	    nim->nifti_type = NIFTI_FTYPE_NIFTI1_1;
	    if(((nim->fname = nifti_strdup(gvnFileName)) == NULL) ||
	       ((nim->iname = nifti_strdup(gvnFileName)) == NULL))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzIVertex2 sz,
			org;
	    sz.vtX = nim->nx;
	    sz.vtY = nim->ny;
	    org.vtX = bBox.xMin;
	    org.vtY = bBox.yMin;
	    errNum = WlzToArray2D(&datAry, obj, sz, org, 0, gType);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    nim->data = *(WlzUByte **)datAry;
	    errno = 0;
	    nifti_image_write(nim);
	    if(errno != 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    nim->data = NULL;
	  }
	  nifti_image_free(nim);
	  (void )Alc2Free(datAry);
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	{
	  WlzIBox3 bBox;
	  void	  ***datAry = NULL; 

	  bBox.xMin = obj->domain.p->kol1;
	  bBox.xMax = obj->domain.p->lastkl;
	  bBox.yMin = obj->domain.p->line1;
	  bBox.yMax = obj->domain.p->lastln;
	  bBox.zMin = obj->domain.p->plane1;
	  bBox.zMax = obj->domain.p->lastpl;
	  if((nim = (nifti_image *)AlcCalloc(1, sizeof(nifti_image))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
            nim->dim[0] = nim->ndim = 3;
	    nim->dim[1] = nim->nx = bBox.xMax - bBox.xMin + 1;
	    nim->dim[2] = nim->ny = bBox.yMax - bBox.yMin + 1;
	    nim->dim[3] = nim->nz = bBox.zMax - bBox.zMin + 1;
	    nim->dim[4] = nim->nt = 1;
	    nim->dim[5] = nim->nu = 1;
	    nim->dim[6] = nim->nv = 1;
	    nim->dim[7] = nim->nw = 1;
	    nim->nvox = nim->nx * nim->ny * nim->nz;
	    nim->datatype = nDType;
	    nim->pixdim[1] = nim->dx = obj->domain.p->voxel_size[0];
	    nim->pixdim[2] = nim->dy = obj->domain.p->voxel_size[1];
	    nim->pixdim[3] = nim->dz = obj->domain.p->voxel_size[2];
	    nim->pixdim[4] = nim->dt = 1.0;
	    nim->pixdim[5] = nim->du = 1.0;
	    nim->pixdim[6] = nim->dv = 1.0;
	    nim->pixdim[7] = nim->dw = 1.0;
	    nim->scl_slope = 1.0;
	    nim->scl_inter = 0.0;
	    nim->sform_code = 1;
	    nim->qto_xyz.m[0][0] = 1.0;
	    nim->qto_xyz.m[1][1] = 1.0;
	    nim->qto_xyz.m[2][2] = 1.0;
	    nim->qto_xyz.m[3][3] = 1.0;
	    nim->qto_ijk.m[0][0] = 1.0;
	    nim->qto_ijk.m[1][1] = 1.0;
	    nim->qto_ijk.m[2][2] = 1.0;
	    nim->qto_ijk.m[3][3] = 1.0;
	    nim->sto_xyz.m[0][0] = 1.0; nim->sto_xyz.m[0][3] = -(bBox.xMin);
	    nim->sto_xyz.m[1][1] = 1.0; nim->sto_xyz.m[1][3] = -(bBox.yMin);
	    nim->sto_xyz.m[2][2] = 1.0; nim->sto_xyz.m[2][3] = -(bBox.zMin);
	    nim->sto_xyz.m[3][3] = 1.0;
	    nim->sto_ijk.m[0][0] = 1.0; nim->sto_ijk.m[0][3] = bBox.xMin;
	    nim->sto_ijk.m[1][1] = 1.0; nim->sto_ijk.m[1][3] = bBox.yMin;
	    nim->sto_ijk.m[2][2] = 1.0; nim->sto_ijk.m[2][3] = bBox.zMin;
	    nim->sto_ijk.m[3][3] = 1.0;
	    nifti_mat44_to_quatern(nim->qto_xyz,
	                           &(nim->quatern_b),
	                           &(nim->quatern_c),
	                           &(nim->quatern_d),
	                           &(nim->qoffset_x),
	                           &(nim->qoffset_y),
	                           &(nim->qoffset_z),
				   NULL, NULL, NULL, &(nim->qfac));
	    nifti_datatype_sizes(nim->datatype,
	                         &(nim->nbyper), &(nim->swapsize) ) ;
            nim->byteorder = nifti_short_order();
	    nim->nifti_type = NIFTI_FTYPE_NIFTI1_1;
	    if(((nim->fname = nifti_strdup(gvnFileName)) == NULL) ||
	       ((nim->iname = nifti_strdup(gvnFileName)) == NULL))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzIVertex3 sz,
			org;
	    sz.vtX = nim->nx;
	    sz.vtY = nim->ny;
	    sz.vtZ = nim->nz;
	    org.vtX = bBox.xMin;
	    org.vtY = bBox.yMin;
	    org.vtZ = bBox.zMin;
	    errNum = WlzToArray3D(&datAry, obj, sz, org, 0, gType);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    nim->data = **(WlzUByte ***)datAry;
	    errno = 0;
	    nifti_image_write(nim);
	    if(errno != 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    nim->data = NULL;
	  }
	  nifti_image_free(nim);
	  (void )Alc3Free(datAry);
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzExtFF
* \brief	Creates a 2D domain object with values from the given NIfTI
* 		image knowing that the NIfTI image is 2D. A compound object
* 		may be created (as determined by WlzEffNiftiToWlzType())
* 		in other cases a 2D domain object is created. The q and s
* 		forms transforms of the NIfTI image are not applied and
* 		neither is the pixel value scaling.
* \param	nim			NIfTI image.
* 					WlzEffReadObjNifti().
* \param	dstErr			Destination error code, may be NULL.
*/
static WlzObject *WlzEffNiftiToObj2D(nifti_image *nim,
				     WlzErrorNum *dstErr)
{
  int		id0,
  		nPPO,
		nVPP;
  size_t	wBPP;
  WlzIVertex2	sz;
  WlzGreyType	wGType;
  WlzObject	*obj = NULL;
  WlzObject	*objs[2];
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = 0;
  objs[0] = objs[1] = NULL;
  wBPP = WlzGreySize(wGType);
  errNum = WlzEffNiftiToWlzType(nim, &nVPP, &wGType);
  if(errNum == WLZ_ERR_NONE)
  {
    sz.vtX = nim->dim[1];
    sz.vtY = nim->dim[2];
    errNum = WlzValueConvertPixel(&bgdV, bgdV, wGType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nPPO = sz.vtX * sz.vtY;
    for(id0 = 0; (errNum == WLZ_ERR_NONE) && (id0 < nVPP); ++id0)
    {
      WlzGreyP	nGP,
      		wGP;

      nGP.v = nim->data;
      if((wGP.v = AlcMalloc(nPPO * wBPP)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	objs[id0] = WlzEffNiftiCopyToObj2D(nGP, nim->datatype, nim->nbyper,
					   wGP, wGType, nVPP, bgdV, sz, id0,
					   &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  AlcFree(wGP.v);
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nVPP > 1)
    {
      (void )WlzAssignObject(objs[0], NULL);
      (void )WlzAssignObject(objs[1], NULL);
      obj = (WlzObject *)WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 2, objs,
                                              WLZ_2D_DOMAINOBJ, &errNum);
    }
    else
    {
      obj = objs[0];
      objs[0] = NULL;
    }
  }
  (void )WlzFreeObj(objs[0]);
  (void )WlzFreeObj(objs[1]);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(obj);
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzExtFF
* \brief	Creates a 2D domain object with values from the given NIfTI
* 		image knowing that the NIfTI image is 2D. A compound object
* 		may be created (as determined by WlzEffNiftiToWlzType())
* 		in other cases a 2D domain object is created. The q and s
* 		forms transforms of the NIfTI image are not applied and
* 		neither is the pixel value scaling.
* \param	nim			NIfTI image.
* \param	dstErr			Destination error code, may be NULL.
*/
static WlzObject *WlzEffNiftiToObj3D(nifti_image *nim,
				     WlzErrorNum *dstErr)
{
  int		id0,
		nVPP;
  WlzIVertex3	sz;
  WlzGreyType	wGType;
  WlzObject	*obj = NULL;
  WlzObject	*objs[2];
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = 0;
  objs[0] = objs[1] = NULL;
  errNum = WlzEffNiftiToWlzType(nim, &nVPP, &wGType);
  if(errNum == WLZ_ERR_NONE)
  {
    sz.vtX = nim->dim[1];
    sz.vtY = nim->dim[2];
    sz.vtZ = nim->dim[3];
    errNum = WlzValueConvertPixel(&bgdV, bgdV, wGType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(id0 = 0; (errNum == WLZ_ERR_NONE) && (id0 < nVPP); ++id0)
    {
      WlzGreyP	nGP;

      nGP.v = nim->data;
      objs[id0] = WlzEffNiftiCopyToObj3D(nGP, nim->datatype, nim->nbyper,
					 wGType, nVPP, bgdV, sz, id0,
					 &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        objs[id0]->domain.p->voxel_size[0] = nim->pixdim[1];
        objs[id0]->domain.p->voxel_size[1] = nim->pixdim[2];
        objs[id0]->domain.p->voxel_size[2] = nim->pixdim[3];
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nVPP > 1)
    {
      (void )WlzAssignObject(objs[0], NULL);
      (void )WlzAssignObject(objs[1], NULL);
      obj = (WlzObject *)WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 2, objs,
                                              WLZ_2D_DOMAINOBJ, &errNum);
    }
    else
    {
      obj = objs[0];
      objs[0] = NULL;
    }
  }
  (void )WlzFreeObj(objs[0]);
  (void )WlzFreeObj(objs[1]);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(obj);
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzExtFF
* \brief	Copies data from the NIfTI image for the given channel
* 		to a new Woolz object.
* \param	nGP			NIfTI data pointer.
* \param	nDT			NIfTI data type.
* \param	nBPP			Bytes per pixel in the NIfTI image.
* \param	wGP			Woolz data pointer.
* \param	wGType			Woolz object grey type.
* \param	nVPP			Number of basic values of wGType per
* 					NIfTI pixel.
* \param	bgdV			Background value for object.
* \param	sz			Size of the object.
* \param	idx			Index of the channel.
* \param	dstErr			Destination error code, may be NULL.
*/
static WlzObject *WlzEffNiftiCopyToObj2D(WlzGreyP nGP, int nDT,
					 int nBPP, WlzGreyP wGP,
					 WlzGreyType wGType, int nVPP,
					 WlzPixelV bgdV, WlzIVertex2 sz,
					 int idx, WlzErrorNum *dstErr)
{
  int		idY;
  size_t	wBPP;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  wBPP = WlzGreySize(wGType);
  obj = WlzMakeRect(0, sz.vtY - 1, 0, sz.vtX - 1,
                    wGType, wGP.inp, bgdV, NULL, NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    for(idY = 0; idY < sz.vtY; ++idY)
    {
      size_t	offY;
      WlzGreyP  nGP2,
      		wGP2;

      offY = idY * sz.vtX;
      nGP2.ubp = (WlzUByte *)(nGP.v) + nBPP * offY;
      wGP2.ubp = wGP.ubp + wBPP * offY;
      WlzEffNiftiCopyToObj1D(wGP2, nGP2, nDT, nVPP, nBPP, sz.vtX, idx);
    }
  }
  if((errNum != WLZ_ERR_NONE) && (obj != NULL))
  {
    obj->values.r->values.v = NULL;
    (void )WlzFreeObj(obj);
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzExtFF
* \brief	Copies data from the NIfTI image for the given channel
* 		to a new Woolz object.
* \param	nGP			NIfTI data pointer.
* \param	nDT			NIfTI data type.
* \param	nBPP			Bytes per pixel in the NIfTI image.
* \param	wGType			Woolz object grey type.
* \param	nVPP			Number of basic values of wGType per
* 					NIfTI pixel.
* \param	bgdV			Background value for object.
* \param	sz			Size of the object.
* \param	idx			Index of the channel.
* \param	dstErr			Destination error code, may be NULL.
*/
static WlzObject *WlzEffNiftiCopyToObj3D(WlzGreyP nGP, int nDT, int nBPP, 
					 WlzGreyType wGType, int nVPP,
					 WlzPixelV bgdV, WlzIVertex3 sz,
					 int idx, WlzErrorNum *dstErr)
{
  int		idZ;
  size_t	wBPP;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  wBPP = WlzGreySize(wGType);
  obj = WlzMakeCuboid(0, sz.vtZ - 1, 0, sz.vtY - 1, 0, sz.vtX - 1,
                      wGType, bgdV, NULL, NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    for(idZ = 0; idZ < sz.vtZ; ++idZ)
    {
      int	idY;
      size_t	offZ;
      WlzGreyP	wGP;

      offZ = idZ * sz.vtX * sz.vtY;
      wGP = (obj->values.vox->values + idZ)->r->values;
      for(idY = 0; idY < sz.vtY; ++idY)
      {
	size_t	  offY;
	WlzGreyP  nGP2,
		  wGP2;

        offY = idY * sz.vtX;
	nGP2.ubp = (WlzUByte *)(nGP.v) + nBPP * (offZ + offY);
	wGP2.ubp = wGP.ubp + wBPP * offY;
	WlzEffNiftiCopyToObj1D(wGP2, nGP2, nDT, nVPP, nBPP, sz.vtX, idx);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \ingroup	WlzExtFF
* \brief	Copies a line of values from a NIfTI image to a Woolz image.
* \param	wGP			Woolz image grey pointer.
* \param	nGP			NIfTI image pointer.
* \param	nDT			NIfTI data type.
* \param	nVPP			Number of basic (Woolz grey) values
* 					per pixel.
* \param	nBPP			Number of bytes per pixel in the
* 					NIfTI image.
* \param	nV			Number of (Woolz) values to set.
* \param	idx			Index of the basic (Woolz grey) value
* 					with each NIfTI pixel.
*/
static void	WlzEffNiftiCopyToObj1D(WlzGreyP wGP, WlzGreyP nGP,
				       int nDT, int nVPP, int nBPP,
				       int nV, int idx)
{
  int		idW;

  switch(nDT)
  {
    case NIFTI_TYPE_UINT8:      /* Unsigned byte, keep as UBYTE. */
      WlzValueCopyUByteToUByte(wGP.ubp, nGP.ubp, nV);
      break;
    case NIFTI_TYPE_INT16:      /* Signed short, keep as short. */
      WlzValueCopyShortToShort(wGP.shp, nGP.shp, nV);
      break;
    case NIFTI_TYPE_INT32:      /* Signed int, keep as int. */
      WlzValueCopyIntToInt(wGP.inp, nGP.inp, nV);
      break;
    case NIFTI_TYPE_FLOAT32: 	/* Float, keep as float. */
      WlzValueCopyFloatToFloat(wGP.flp, nGP.flp, nV);
      break;
    case NIFTI_TYPE_COMPLEX64: 	/* Complex float, convert to 2 x float
				   objects. */
      nGP.flp += idx;
      for(idW = 0; idW < nV; ++idW)
      {
        *(wGP.flp + idW) = *(nGP.flp + (nVPP * idW));
      }
      break;
    case NIFTI_TYPE_FLOAT64:	/* Double, keep as double. */
      WlzValueCopyDoubleToDouble(wGP.dbp, nGP.dbp, nV);
      break;
    case NIFTI_TYPE_RGB24:	/* RGB, promote to RGBA. */
      for(idW = 0; idW < nV; ++idW)
      {
	WlzUInt	r;

        WLZ_RGBA_RED_SET(r,   *(nGP.ubp + 0));
        WLZ_RGBA_GREEN_SET(r, *(nGP.ubp + 1));
        WLZ_RGBA_BLUE_SET(r,  *(nGP.ubp + 2));
        WLZ_RGBA_ALPHA_SET(r, 255);
	nGP.ubp += nBPP;
      }
      break;
    case NIFTI_TYPE_INT8:	/* Signed byte, promote to short. */
      for(idW = 0; idW < nV; ++idW)
      {
        *(wGP.shp + idW) = *(char *)(nGP.ubp + idW);
      }
      break;
    case NIFTI_TYPE_UINT16:	/* Unsigned short, promote to int. */
      for(idW = 0; idW < nV; ++idW)
      {
        *(wGP.inp + idW) = *(unsigned short *)(nGP.shp + idW);
      }
      break;
    case NIFTI_TYPE_UINT32:	/* Unsigned int, promote to long long. */
      for(idW = 0; idW < nV; ++idW)
      {
        *(wGP.lnp + idW) = *(unsigned int *)(nGP.shp + idW);
      }
      break;
    case NIFTI_TYPE_INT64:	/* Long long, keep long long. */
      for(idW = 0; idW < nV; ++idW)
      {
        *(wGP.lnp + idW) = *(nGP.lnp + idW);
      }
      break;
    case NIFTI_TYPE_UINT64:	/* Unsigned long long, treat as long long. */
      for(idW = 0; idW < nV; ++idW)
      {
        *(wGP.lnp + idW) = WLZ_CLAMP(*(nGP.lnp + idW), LLONG_MIN, LLONG_MAX);
      }
      break;
    case NIFTI_TYPE_FLOAT128:   /* Long double, error. */
      break;
    case NIFTI_TYPE_COMPLEX128:	/* Complex double, convert to 2 x double
                                       objects. */
      nGP.dbp += idx;
      for(idW = 0; idW < nV; ++idW)
      {
        *(wGP.dbp + idW) = *(nGP.dbp + (nVPP * idW));
      }
      break;
    case NIFTI_TYPE_COMPLEX256:	/* Complex long double, error. */
      break;
    case NIFTI_TYPE_RGBA32:	/* RGBA, keep as RGBA. */
      WlzValueCopyRGBAToRGBA(wGP.rgbp, nGP.rgbp, nV);
      break;
    default:
      break;
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Computes Woolz pixel types and object count from the NIfTI
*		data type. Some NIfTI images need to be represented
*		by a pair of Woolz objects, eg images with complex image
*		values.
* \param	nim			NIfTI image.
* \param	dstNVPP		Destination pointer for the number of
* 					Woolz grey values per NIfTI pixel,
* 					must not be NULL.
* \param	dstWGType		Destination pointer for the Woolz
* 					grey type, must not be NULL.
*/
static WlzErrorNum WlzEffNiftiToWlzType(nifti_image *nim, int *dstNVPP,
				        WlzGreyType *dstWGType)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(nim->datatype)
  {
    case NIFTI_TYPE_UINT8:	/* Unsigned byte, keep UBYTE. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_UBYTE;
      break;
    case NIFTI_TYPE_INT16:      /* Signed short, keep as short. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_SHORT;
      break;
    case NIFTI_TYPE_INT32:      /* Signed int, keep as int. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_INT;
      break;
    case NIFTI_TYPE_FLOAT32: 	/* Float, keep as float. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_FLOAT;
      break;
    case NIFTI_TYPE_COMPLEX64: 	/* Complex float, convert to 2 x float
				   objects. */
      *dstNVPP = 2;
      *dstWGType = WLZ_GREY_FLOAT;
      break;
    case NIFTI_TYPE_FLOAT64:	/* Double, keep as double. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_DOUBLE;
      break;
    case NIFTI_TYPE_RGB24:	/* RGB, promote to RGBA. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_RGBA;
      break;
    case NIFTI_TYPE_INT8:      /* Signed byte, promote to short. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_SHORT;
      break;
    case NIFTI_TYPE_UINT16:	/* Unsigned short, promote to int. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_SHORT;
      break;
    case NIFTI_TYPE_UINT32:	/* Unsigned int, promote to long long. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_LONG;
      break;
    case NIFTI_TYPE_INT64:	/* Long long, keep long long. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_LONG;
      break;
    case NIFTI_TYPE_UINT64:	/* Unsigned long long, treat as long long. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_LONG;
      break;
    case NIFTI_TYPE_FLOAT128:   /* Long double, error. */
      errNum = WLZ_ERR_GREY_TYPE;
      break;
    case NIFTI_TYPE_COMPLEX128:	/* Complex double, convert to 2 x double
                                       objects. */
      *dstNVPP = 2;
      *dstWGType = WLZ_GREY_DOUBLE;
      break;
    case NIFTI_TYPE_COMPLEX256:	/* Complex long double, error. */
      errNum = WLZ_ERR_GREY_TYPE;
      break;
    case NIFTI_TYPE_RGBA32:	/* RGBA, keep as RGBA. */
      *dstNVPP = 1;
      *dstWGType = WLZ_GREY_RGBA;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return	Formated ASCII string.
* \ingroup	WlzExtFF
* \brief	Allocates and fills an ASCII string with the NIfTI image
* 		information and description.
* \param	nim			NIfTI image.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static char 	*WlzEffNiftiInfoToStr(nifti_image *nim, WlzErrorNum *dstErr)
{
  char		*s;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	maxSLen = 4096;

  if((s = (char *)AlcMalloc(sizeof(char) * maxSLen)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    (void )sprintf(s,
                   "dim=%d,%d,%d,%d,%d,%d,%d,%d\n"
		   "intent_code=%d\n"
		   "intent_name=\"%s\"\n"
		   "intent_p1=%g\n"
		   "intent_p2=%g\n"
		   "intent_p3=%g\n"
		   "pixdim=%g,%g,%g,%g,%g,%g,%g,%g\n"
		   "qform_code=%d\n"
		   "qto_ijk=%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n"
		   "sform_code=%d\n"
		   "slice_code=%d\n"
		   "slice_start=%d\n"
		   "slice_end=%d\n"
		   "slice_duration=%g\n"
		   "cal_min=%g\n"
		   "cal_max=%g\n"
		   "toffset=%g\n"
		   "Description = \"%s\"\n",
		   nim->dim[0], nim->dim[1], nim->dim[2], nim->dim[3],
		   nim->dim[4], nim->dim[5], nim->dim[6], nim->dim[7],
		   nim->intent_code, nim->intent_name,
		   nim->intent_p1, nim->intent_p2, nim->intent_p3,
		   nim->pixdim[0], nim->pixdim[1],
		   nim->pixdim[2], nim->pixdim[3],
		   nim->pixdim[4], nim->pixdim[5],
		   nim->pixdim[6], nim->pixdim[7],
		   nim->qform_code,
		   nim->qto_ijk.m[0][0], nim->qto_ijk.m[0][1],
		   nim->qto_ijk.m[0][2], nim->qto_ijk.m[0][3],
		   nim->qto_ijk.m[1][0], nim->qto_ijk.m[1][1],
		   nim->qto_ijk.m[1][2], nim->qto_ijk.m[1][3],
		   nim->qto_ijk.m[2][0], nim->qto_ijk.m[2][1],
		   nim->qto_ijk.m[2][2], nim->qto_ijk.m[2][3],
		   nim->qto_ijk.m[3][0], nim->qto_ijk.m[3][1],
		   nim->qto_ijk.m[3][2], nim->qto_ijk.m[3][3],
		   nim->sform_code,
		   nim->slice_code,
		   nim->slice_start,
		   nim->slice_end,
		   nim->slice_duration,
		   nim->cal_min,
		   nim->cal_max,
		   nim->toffset,
		   nim->descrip);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(s);
}

/*!
* \return	NIfTI data type code.
* \ingroup	WlzExtFF
* \brief	Converts the given Woolz grey type to the matching NIfTI data
* 		type.
* \param	wGType			Woolz grey type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static int 	WlzEffNiftiFromWlzGType(WlzGreyType wGType, WlzErrorNum *dstErr)
{
  int 		nDType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(wGType)
  {
    case WLZ_GREY_UBYTE:
      nDType = NIFTI_TYPE_UINT8;
      break;
    case WLZ_GREY_SHORT:
      nDType = NIFTI_TYPE_INT16;
      break;
    case WLZ_GREY_INT:
      nDType = NIFTI_TYPE_INT32;
      break;
    case WLZ_GREY_RGBA:
      nDType = DT_RGBA32;
      break;
    case WLZ_GREY_FLOAT:
      nDType = DT_FLOAT32;
      break;
    case WLZ_GREY_DOUBLE:
      nDType = DT_FLOAT64;
      break;
    default:
      nDType = DT_NONE;
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nDType);
}
