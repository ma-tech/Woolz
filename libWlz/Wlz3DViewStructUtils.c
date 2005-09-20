#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/Wlz3DViewStructUtils.c
* \author       Richard Baldock, Bill Hill
* \date         October 2001
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
* \brief	Utility functions associated with 3D views.
* \ingroup	WlzSectionTransform
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static double viewStructAtan2(
  double x,
  double y)
{
  if( (fabs(y) < 1.e-6) && (fabs(x) < 1.e-6) ){
    return 0.0;
  }
  return atan2(x, y);
}

/*!
* \return				The 3D view or NULL on error.
* \ingroup      WlzSectionTransform
* \brief	Reads a 3D section structure from a file.
* \param	fp			Input file pointer.
* \param	dstErr			Destination pointer for an error
*					code, may be NULL.
*/
WlzThreeDViewStruct *WlzRead3DViewStruct(
  FILE		*fp,
  WlzErrorNum	*dstErr)
{
  if( dstErr ){
    *dstErr = WLZ_ERR_UNSPECIFIED;
  }
  return NULL;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Writes a 3D section structure to a file.
* \param	fp			Output file pointer.
* \param	viewstr			Given view.
*/
WlzErrorNum WlzWrite3DViewStruct(
  FILE			*fp,
  WlzThreeDViewStruct	*viewstr)
{
  return WLZ_ERR_UNSPECIFIED ;
}

/*!
* \return				A new view.
* \ingroup      WlzSectionTransform
* \brief	Allocates and intialises a 3D view. The transform
*		look up tables are not allocated.
* \param	type			Only WLZ_3D_VIEW_STRUCT is
*					currently allowed.
* \param	dstErr			Destination pointer for an error
*					code, may be NULL.
*/
WlzThreeDViewStruct *WlzMake3DViewStruct(
  WlzObjectType	type,
  WlzErrorNum	*dstErr)
{
  WlzThreeDViewStruct	*viewStr;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzAffineTransform	*trans = NULL;

  /* check type - in principle that is, when defined!! */
  switch( type ){

  case WLZ_3D_VIEW_STRUCT:
    break;

  default:
    errNum = WLZ_ERR_OBJECT_TYPE;
    viewStr = NULL;

  }

  /* allocate space */
  if( errNum == WLZ_ERR_NONE ){
    if( (viewStr = (WlzThreeDViewStruct *)
	 AlcCalloc(1, sizeof(WlzThreeDViewStruct))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  /* set default values */
  if( errNum == WLZ_ERR_NONE ){
    viewStr->type 	= type;
    viewStr->linkcount 	= 0;
    viewStr->freeptr 	= NULL;
    viewStr->initialised= 0;
    viewStr->fixed.vtX 	= 0.0;
    viewStr->fixed.vtY 	= 0.0;
    viewStr->fixed.vtZ 	= 0.0;
    viewStr->theta 	= 0.0;
    viewStr->phi 	= 0.0;
    viewStr->zeta 	= 0.0;
    viewStr->dist 	= 0.0;
    viewStr->scale 	= 1.0;
    viewStr->voxelSize[0] = 1.0;
    viewStr->voxelSize[1] = 1.0;
    viewStr->voxelSize[2] = 1.0;
    viewStr->voxelRescaleFlg = 0;
    viewStr->interp = WLZ_INTERPOLATION_NEAREST;
    viewStr->view_mode	= WLZ_STATUE_MODE;
    viewStr->up.vtX 	= 0.0;
    viewStr->up.vtY 	= 0.0;
    viewStr->up.vtZ 	= 1.0;
    viewStr->ref_obj 	= NULL;
    /*AlcDouble2Malloc(&viewStr->rotation, 3, 3);*/
    trans = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, &errNum);
    viewStr->trans = WlzAssignAffineTransform(trans, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return viewStr ;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Frees a view.
* \param	viewStr			Given view to free.
*/
WlzErrorNum WlzFree3DViewStruct(
  WlzThreeDViewStruct *viewStr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check pointer */
  if( viewStr == NULL ){
    return WLZ_ERR_DOMAIN_NULL;
  }

  /* check type */

  /* check linkcount */
  if( WlzUnlink(&(viewStr->linkcount), &errNum) )
  {
    /* free luts and rotation matrix if set */
    if( viewStr->initialised && viewStr->freeptr != NULL ){
      AlcFreeStackFree(viewStr->freeptr);
    }
    if( viewStr->ref_obj != NULL ){
      WlzFreeObj( viewStr->ref_obj );
    }

    /* free the affine transform and structure */
    /*AlcDouble2Free(viewStr->rotation);*/
    (void) WlzFreeAffineTransform(viewStr->trans);
    AlcFree((void *) viewStr);
  }

  return errNum;
}

/*!
* \return
* \ingroup      WlzSectionTransform
* \brief	Sets up a rotation matrix from the given Euler
* 		angles.
* \param	rotation		3x3 rotation matrix.
* \param	xsi			Rotation about the z-axis.
* \param	eta			Rotation about the new y-axis.
* \param	zeta			Rotation about the new z-axis.
*/
static void setupEulerRotationMatrix(
  double	**rotation,
  double	xsi,
  double	eta,
  double	zeta)
{
  double cos_x = cos(xsi), sin_x = sin(xsi);
  double cos_e = cos(eta), sin_e = sin(eta);
  double cos_z = cos(zeta), sin_z = sin(zeta);

  rotation[0][0] =  cos_z * cos_e * cos_x - sin_z * sin_x;
  rotation[0][1] =  cos_z * cos_e * sin_x + sin_z * cos_x;
  rotation[1][0] = -sin_z * cos_e * cos_x - cos_z * sin_x;
  rotation[1][1] = -sin_z * cos_e * sin_x + cos_z * cos_x;
  rotation[0][2] = -cos_z * sin_e;
  rotation[1][2] =  sin_z * sin_e;
  rotation[2][0] =  sin_e * cos_x;
  rotation[2][1] =  sin_e * sin_x;
  rotation[2][2] =  cos_e;

  return;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets up the rotation matrix of the given view.
* \param	viewStr			Given view.
*/
static WlzErrorNum setup3DSectionAffineTransform(
  WlzThreeDViewStruct	*viewStr)
{
  double xsi, eta, zeta;
  double cos_t, sin_t, cos_p, sin_p;
  double upp_x, upp_y;
  double fx, fy, fz;
  double **rotation;
  WlzAffineTransform	*tmpTrans, *tr, *tf, *ts;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int i, j;

  /* check pointer */
  if( viewStr == NULL ){
    return WLZ_ERR_OBJECT_NULL;
  }

  /* switch on the mode */
  xsi = viewStr->theta;
  eta = viewStr->phi;
  switch( viewStr->view_mode ){
  case WLZ_STATUE_MODE:
    zeta = -xsi;
    break;

  case WLZ_UP_IS_UP_MODE:
    cos_t = cos(viewStr->theta);
    sin_t = sin(viewStr->theta);
    cos_p = cos(viewStr->phi);
    sin_p = sin(viewStr->phi);
    if( viewStr->voxelRescaleFlg ){
      upp_x = (viewStr->up.vtX*viewStr->voxelSize[0]*cos_p*cos_t +
	       viewStr->up.vtY*viewStr->voxelSize[1]*cos_p*sin_t -
	       viewStr->up.vtZ*viewStr->voxelSize[2]*sin_p);
      upp_y = (-viewStr->up.vtX*viewStr->voxelSize[0]*sin_t +
	       viewStr->up.vtY*viewStr->voxelSize[1]*cos_t);
    }
    else {
      upp_x = (viewStr->up.vtX*cos_p*cos_t + viewStr->up.vtY*cos_p*sin_t -
	       viewStr->up.vtZ*sin_p);
      upp_y = (-viewStr->up.vtX*sin_t + viewStr->up.vtY*cos_t);
    }
    zeta = viewStructAtan2(upp_x, upp_y);
    break;

  case WLZ_ZERO_ZETA_MODE:
    zeta = 0.0;
    break;

  case WLZ_ZETA_MODE:
    zeta = viewStr->zeta;
    break;

  case WLZ_FIXED_LINE_MODE:
    /* this will need fixing -- done */
    cos_t = cos(viewStr->theta);
    sin_t = sin(viewStr->theta);
    cos_p = cos(viewStr->phi);
    sin_p = sin(viewStr->phi);
    fx = viewStr->fixed_2.vtX - viewStr->fixed.vtX;
    fy = viewStr->fixed_2.vtY - viewStr->fixed.vtY;
    fz = viewStr->fixed_2.vtZ - viewStr->fixed.vtZ;
    if( viewStr->voxelRescaleFlg ){
      fx *= viewStr->voxelSize[0];
      fy *= viewStr->voxelSize[1];
      fz *= viewStr->voxelSize[2];
    }
    upp_x = (fx*cos_p*cos_t + fy*cos_p*sin_t - fz*sin_p);
    upp_y = (-fx*sin_t + fy*cos_t);
    zeta = -(viewStr->fixed_line_angle - viewStructAtan2(upp_y, upp_x));
    break;
    
  default:
    return WLZ_ERR_INT_DATA;
  }

  /* get zeta into the range 0 to 2*pi */
  while( zeta > 2 * WLZ_M_PI ){ zeta -= 2 * WLZ_M_PI;}
  while( zeta < 0.0 ){ zeta += 2 * WLZ_M_PI;}
  viewStr->zeta = zeta;

  /* allocate workspace for the augmented transform matrix */
  AlcDouble2Malloc(&rotation, 4, 4);
  setupEulerRotationMatrix(rotation, xsi, eta, zeta);

  /* build the transform matrix */
  rotation[3][3] = 1.0;
  for(i=0; i < 3; i++){
    rotation[i][3] = 0.0;
    rotation[3][i] = 0.0;
  }
  tr = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_3D_AFFINE,
				    rotation, &errNum);

  tf = WlzAffineTransformFromTranslation(WLZ_TRANSFORM_3D_AFFINE,
					 -viewStr->fixed.vtX,
					 -viewStr->fixed.vtY,
					 -viewStr->fixed.vtZ,
					 &errNum);

  /* check for rescaling */
  if( viewStr->voxelRescaleFlg ){
    /* define re-scaling affine transform */
    for(i=0; i < 4; i++){
      for(j=0; j < 4; j++){
	rotation[i][j] = 0.0;
      }
    }
    for(i=0; i < 3; i++){
      rotation[i][i] = viewStr->voxelSize[i];
    }
    rotation[3][3] = 1.0;
    ts = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_3D_AFFINE,
				      rotation, &errNum);
  }
  else {
    ts = WlzAffineTransformFromTranslation(WLZ_TRANSFORM_3D_AFFINE,
					   0.0, 0.0, 0.0, &errNum);
  }
  tmpTrans = WlzAffineTransformProduct(tf, ts, &errNum);
  WlzFreeAffineTransform(tf);
  WlzFreeAffineTransform(ts);
  ts = tmpTrans;
  tmpTrans = WlzAffineTransformProduct(ts, tr, &errNum);
  WlzFreeAffineTransform(ts);
  WlzFreeAffineTransform(tr);

  /* assign to the transform */
  WlzAffineTransformMatrixSet(viewStr->trans, tmpTrans->mat);
  WlzFreeAffineTransform(tmpTrans);
  AlcDouble2Free(rotation);

  return WLZ_ERR_NONE;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets up the transformation look up tables of the given view.
* \param	viewStr			Given view.
*/
static void setupTransformLuts(
  WlzThreeDViewStruct	*viewStr)
{
  int		xp, yp, i;
  WlzAffineTransform	*trans;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzDVertex3	vtx, dst;

  /* get the inverse transform */
  trans = WlzAffineTransformInverse(viewStr->trans, &errNum);

  for(xp= WLZ_NINT(viewStr->minvals.vtX), i=0;
      xp <= WLZ_NINT(viewStr->maxvals.vtX); xp++, i++){
    vtx.vtX = xp;
    vtx.vtY = 0.0;
    vtx.vtZ = viewStr->dist;
    dst = WlzAffineTransformVertexD3(trans, vtx, &errNum);
    viewStr->xp_to_x[i] = dst.vtX;
    viewStr->xp_to_y[i] = dst.vtY;
    viewStr->xp_to_z[i] = dst.vtZ;
  }

  for(yp= WLZ_NINT(viewStr->minvals.vtY), i=0;
      yp <= WLZ_NINT(viewStr->maxvals.vtY); yp++, i++){
    viewStr->yp_to_x[i] = trans->mat[0][1]*yp;
    viewStr->yp_to_y[i] = trans->mat[1][1]*yp;
    viewStr->yp_to_z[i] = trans->mat[2][1]*yp;
  }

  WlzFreeAffineTransform(trans);
  return;
}


#define CHECK_MIN_MAX_VTX(X,Y,Z) \
vtx.vtX = X; vtx.vtY = Y; vtx.vtZ = Z; \
Wlz3DSectionTransformVtx( &vtx, viewStr ); \
if( vtx.vtX < viewStr->minvals.vtX )viewStr->minvals.vtX = vtx.vtX; \
if( vtx.vtX > viewStr->maxvals.vtX )viewStr->maxvals.vtX = vtx.vtX; \
if( vtx.vtY < viewStr->minvals.vtY )viewStr->minvals.vtY = vtx.vtY; \
if( vtx.vtY > viewStr->maxvals.vtY )viewStr->maxvals.vtY = vtx.vtY; \
if( vtx.vtZ < viewStr->minvals.vtZ )viewStr->minvals.vtZ = vtx.vtZ; \
if( vtx.vtZ > viewStr->maxvals.vtZ )viewStr->maxvals.vtZ = vtx.vtZ;

#define CHECK_3D_BB_VTX(X,Y) \
vtx.vtX = X; \
vtx.vtY = Y; \
vtx.vtZ = viewStr->dist; \
Wlz3DSectionTransformInvVtx( &vtx, viewStr ); \
i = WLZ_NINT(vtx.vtX); \
tmpPlanedmn.kol1 = WLZ_MIN(tmpPlanedmn.kol1, i); \
tmpPlanedmn.lastkl = WLZ_MAX(tmpPlanedmn.lastkl, i); \
i = WLZ_NINT(vtx.vtY); \
tmpPlanedmn.line1 = WLZ_MIN(tmpPlanedmn.line1, i); \
tmpPlanedmn.lastln = WLZ_MAX(tmpPlanedmn.lastln, i); \
i = WLZ_NINT(vtx.vtZ); \
tmpPlanedmn.plane1 = WLZ_MIN(tmpPlanedmn.plane1, i); \
tmpPlanedmn.lastpl = WLZ_MAX(tmpPlanedmn.lastpl, i);



/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Initialises a 3D view with respect to the given view parameters
*		and a 3D object.
*
*		Initialises a 3D view with respect to the given view parameters
*		and a 3D object. Initialisation involves calculating the bounds
*		of the section, setting up the look up tables to calculate the
*		section image, setting up the rotation matrix and attaching the
*		object as the reference object for this view. This view
*		structure can be reused for simple changes of the view
*		parameter "dist" but otherwise must be re-initialised.
* \param	viewStr			View to be intialised.
* \param	obj			The 3D object to be sectioned.
*/
WlzErrorNum WlzInit3DViewStruct(
  WlzThreeDViewStruct	*viewStr,
  WlzObject		*obj)
{
  double		*tDP0;
  WlzPlaneDomain	*planedmn, tmpPlanedmn;
  WlzDVertex3		vtx;
  unsigned int		widthp, heightp;
  int			i;
  WlzErrorNum		dstErr;
  

  /* check the pointers and types */
  if( viewStr == NULL || obj == NULL ){
    return WLZ_ERR_OBJECT_NULL;
  }
  switch( obj->type ){

  case WLZ_3D_DOMAINOBJ: /* 3D object to be sectioned */
    planedmn = obj->domain.p;
    if( planedmn == NULL ){
      return WLZ_ERR_DOMAIN_NULL;
    }
    if( planedmn->type != WLZ_PLANEDOMAIN_DOMAIN ){
      return WLZ_ERR_DOMAIN_TYPE;
    }

    /* put in the voxel sizes, rescaled to min 1.0 */
    viewStr->voxelSize[0] = planedmn->voxel_size[0];
    viewStr->voxelSize[1] = planedmn->voxel_size[1];
    viewStr->voxelSize[2] = planedmn->voxel_size[2];

    /* set the rotation matrix */
    (void) setup3DSectionAffineTransform(viewStr);
    break;

  case WLZ_2D_DOMAINOBJ:	/* assume this is the section */
    if( obj->domain.i == NULL ){
      return WLZ_ERR_DOMAIN_TYPE;
    }
    /* get the 3D bounding box */
    (void) setup3DSectionAffineTransform(viewStr);
    vtx.vtX = obj->domain.i->kol1;
    vtx.vtY = obj->domain.i->line1;
    vtx.vtZ = viewStr->dist;
    Wlz3DSectionTransformInvVtx( &vtx, viewStr );
    tmpPlanedmn.kol1 = WLZ_NINT(vtx.vtX);
    tmpPlanedmn.line1 = WLZ_NINT(vtx.vtY);
    tmpPlanedmn.plane1 = WLZ_NINT(vtx.vtZ);
    tmpPlanedmn.lastkl = tmpPlanedmn.kol1;
    tmpPlanedmn.lastln = tmpPlanedmn.line1;
    tmpPlanedmn.lastpl = tmpPlanedmn.plane1;
    CHECK_3D_BB_VTX(obj->domain.i->kol1, obj->domain.i->lastln);
    CHECK_3D_BB_VTX(obj->domain.i->lastkl, obj->domain.i->line1);
    CHECK_3D_BB_VTX(obj->domain.i->lastkl, obj->domain.i->lastln);
    tmpPlanedmn.kol1 -= 1;
    tmpPlanedmn.line1 -= 1;
    tmpPlanedmn.plane1 -= 1;
    tmpPlanedmn.lastkl += 1;
    tmpPlanedmn.lastln += 1;
    tmpPlanedmn.lastpl += 1;
    planedmn = &tmpPlanedmn;

    break;

  default:
    return WLZ_ERR_OBJECT_TYPE;

  }

  /* if intialised then free existing luts */
  if( viewStr->initialised ){
    if( viewStr->freeptr ){
      AlcFreeStackFree( viewStr->freeptr );
    }
  }

  /* transform one corner of the bounding box to initialise min and max
     coordinate values */
  vtx.vtX = planedmn->kol1;
  vtx.vtY = planedmn->line1;
  vtx.vtZ = planedmn->plane1;
  Wlz3DSectionTransformVtx( &vtx, viewStr );
  viewStr->minvals.vtX = viewStr->maxvals.vtX = vtx.vtX;
  viewStr->minvals.vtY = viewStr->maxvals.vtY = vtx.vtY;
  viewStr->minvals.vtZ = viewStr->maxvals.vtZ = vtx.vtZ;
  CHECK_MIN_MAX_VTX(planedmn->kol1, planedmn->line1, planedmn->lastpl);
  CHECK_MIN_MAX_VTX(planedmn->kol1, planedmn->lastln, planedmn->plane1);
  CHECK_MIN_MAX_VTX(planedmn->kol1, planedmn->lastln, planedmn->lastpl);
  CHECK_MIN_MAX_VTX(planedmn->lastkl, planedmn->line1, planedmn->plane1);
  CHECK_MIN_MAX_VTX(planedmn->lastkl, planedmn->line1, planedmn->lastpl);
  CHECK_MIN_MAX_VTX(planedmn->lastkl, planedmn->lastln, planedmn->plane1);
  CHECK_MIN_MAX_VTX(planedmn->lastkl, planedmn->lastln, planedmn->lastpl);

  /* find range of values required, and allocate space for the LUT's */
  widthp  = WLZ_NINT(viewStr->maxvals.vtX) -
    WLZ_NINT(viewStr->minvals.vtX) + 1;
  heightp = WLZ_NINT(viewStr->maxvals.vtY) -
    WLZ_NINT(viewStr->minvals.vtY) + 1;
  AlcDouble1Malloc(&tDP0, 3*widthp + 3*heightp);
  viewStr->freeptr = AlcFreeStackPush(NULL, tDP0, NULL);
  viewStr->xp_to_x = tDP0;
  viewStr->xp_to_y = viewStr->xp_to_x + widthp;
  viewStr->xp_to_z = viewStr->xp_to_y + widthp;
  viewStr->yp_to_x = viewStr->xp_to_z + widthp;
  viewStr->yp_to_y = viewStr->yp_to_x + heightp;
  viewStr->yp_to_z = viewStr->yp_to_y + heightp;

  /* set LUT values using current parameters */
  setupTransformLuts( viewStr );

  /* set initialised and attach the reference object */
  viewStr->initialised = 1;
  if( viewStr->ref_obj != obj ){
    if( viewStr->ref_obj ){
      WlzFreeObj( viewStr->ref_obj );
    }
    viewStr->ref_obj = WlzAssignObject( obj, &dstErr );
  }

  return WLZ_ERR_NONE;
}

/*!
* \return
* \ingroup      WlzSectionTransform
* \brief	Transforms a 3D vertex using the section transform
*		overwriting the vertex values.
* \param	vtx			Given vertex, values are overwritten.
* \param	viewStr			Given view.
*/
WlzErrorNum Wlz3DSectionTransformVtx(
  WlzDVertex3		*vtx,
  WlzThreeDViewStruct	*viewStr)
{
    WlzDVertex3 dst;
    WlzErrorNum	errNum=WLZ_ERR_NONE;

    dst = WlzAffineTransformVertexD3(viewStr->trans, *vtx, &errNum);
    *vtx = dst;

    return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Transforms a 3D vertex using the section transform.
* \param	viewStr			Given view.
* \param	vtx			Vertex to be transformed.
* \param	dstVtx			destination pointer for the transformed
* 					vertex.
*/
WlzErrorNum Wlz3DSectionTransformVtxR(
  WlzThreeDViewStruct	*viewStr,
  WlzDVertex3		vtx,
  WlzDVertex3		*dstVtx)
{
  WlzDVertex3	new;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if (dstVtx == NULL){
	  errNum = WLZ_ERR_PARAM_NULL;
  }
  else if( viewStr == NULL ){
	  errNum = WLZ_ERR_PARAM_NULL;
  }
  else {
    new.vtX = vtx.vtX;
    new.vtY = vtx.vtY;
    new.vtZ = vtx.vtZ;
    errNum = Wlz3DSectionTransformVtx(&new, viewStr);
    dstVtx->vtX = new.vtX;
    dstVtx->vtY = new.vtY;
    dstVtx->vtZ = new.vtZ;
  }
  
  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Inverse transforms a 3D vertex using the section transform,
*		overwriting the vertex values.
* \param	vtx			Given 3D vertex, values are
* 					overwritten.
* \param	viewStr			Given view.
*/
WlzErrorNum Wlz3DSectionTransformInvVtx(
  WlzDVertex3		*vtx,
  WlzThreeDViewStruct	*viewStr)
{
    WlzDVertex3		dst;
    WlzAffineTransform	*trans;
    WlzErrorNum		errNum=WLZ_ERR_NONE;

    if( trans = WlzAffineTransformInverse(viewStr->trans, &errNum) ){
      dst = WlzAffineTransformVertexD3(trans, *vtx, &errNum);
      *vtx = dst;
      WlzFreeAffineTransform(trans);
    }

    return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Inverse transforms a 3D vertex using the section transform.
* \param	viewStr			Given view.
* \param	vtx			Vertex to be transformed.
* \param	dstVtx			Destination pointer for the transformed
* 					vertex.
*/
WlzErrorNum Wlz3DSectionTransformInvVtxR(
  WlzThreeDViewStruct	*viewStr,
  WlzDVertex3		vtx,
  WlzDVertex3		*dstVtx)
{
  WlzDVertex3	new;
  WlzErrorNum	errNum;

  if (dstVtx == NULL){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if( viewStr == NULL ){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else {
    new.vtX = vtx.vtX;
    new.vtY = vtx.vtY;
    new.vtZ = vtx.vtZ;
    errNum = Wlz3DSectionTransformInvVtx(&new, viewStr);
    dstVtx->vtX = new.vtX;
    dstVtx->vtY = new.vtY;
    dstVtx->vtZ = new.vtZ;
  }
  
  return errNum;
}

/*!
* \return				Always returns WLZ_ERR_NONE because
*					this fuction does no error checking!
* \ingroup      WlzSectionTransform
* \brief	Increments the distance parameter of a 3D view resetting the
* 		look up tables as required. This is provided because changing
* 		the distance does not require the rotation matrix and all the
* 		look up tables only a single multiply and add per look up
* 		table entry instead of 2 multiplies and three adds plus all the
* 		y look up tables and the rotation matrix with many trigonometry
*		calculations.
* \param	viewStr			Given view.
* \param	incr			The increment distance.
*/
WlzErrorNum Wlz3DSectionIncrementDistance(
  WlzThreeDViewStruct	*viewStr,
  double		incr)
{
  int		xp;
  WlzAffineTransform	*trans;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* get the inverse transform */
  if( trans = WlzAffineTransformInverse(viewStr->trans, &errNum) ){
    for(xp=0; xp <=  WLZ_NINT(viewStr->maxvals.vtX) -
	  WLZ_NINT(viewStr->minvals.vtX); xp++){
      viewStr->xp_to_x[xp] += trans->mat[0][2]*incr;
      viewStr->xp_to_y[xp] += trans->mat[1][2]*incr;
      viewStr->xp_to_z[xp] += trans->mat[2][2]*incr;
    }
    errNum = WlzFreeAffineTransform(trans);
    viewStr->dist += incr;
  }

  return errNum;
}

/*!
* \return				The point of intersection.
* \ingroup      WlzSectionTransform
* \brief	Finds a point on the line of intersection of two 3D views. The
* 		point is returned in the coordinate system of the first view.
*		This function together with Wlz3DViewGetIntersectionAngle() is
*		useful for finding the line of intersection within one of the
*		planes.
* \param	v1			First view.
* \param	v2			Second view.
* \param	dstErr			Destination pointer for an error
*					code, may be NULL.
*/
WlzDVertex2 Wlz3DViewGetIntersectionPoint(
  WlzThreeDViewStruct *v1,
  WlzThreeDViewStruct *v2,
  WlzErrorNum	*dstErr)
{
  WlzDVertex2	rtnVtx;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzAffineTransform	*t1, *t2;
  double	*a[3], ap[3][3], b[3];
  double	dp1=0.0, dp12=0.0;
  int		i;
  double	d1 = v1->dist;
  double	d2 = v2->dist;

  /* check view structs - non-NULL and initialised */
  if( (v1 == NULL) || (v2 == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if( !v1->initialised || !v2->initialised ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else {
    if( t1 = WlzAffineTransformInverse(v1->trans, &errNum) ){
      if( t2 = WlzAffineTransformInverse(v2->trans, &errNum) ){
	for(i=0; i < 3; i++){
	  dp1 += t1->mat[i][2] * t1->mat[i][2];
	  dp12 += t1->mat[i][2] * t2->mat[i][2];
	}
	ap[0][0] = t2->mat[0][0];
	ap[1][0] = t2->mat[1][0];
	ap[2][0] = t2->mat[2][0];
	ap[0][1] = t2->mat[0][1];
	ap[1][1] = t2->mat[1][1];
	ap[2][1] = t2->mat[2][1];
	ap[0][2] = -(dp12*t1->mat[0][2] - dp1*t2->mat[0][2]);
	ap[1][2] = -(dp12*t1->mat[1][2] - dp1*t2->mat[1][2]);
	ap[2][2] = -(dp12*t1->mat[2][2] - dp1*t2->mat[2][2]);
	a[0] = &ap[0][0];
	a[1] = &ap[1][0];
	a[2] = &ap[2][0];

	b[0] = t1->mat[0][2]*d1 + t1->mat[0][3]
	  - t2->mat[0][2]*d2 - t2->mat[0][3];
	b[1] = t1->mat[1][2]*d1 + t1->mat[1][3]
	  - t2->mat[1][2]*d2 - t2->mat[1][3];
	b[2] = t1->mat[2][2]*d1 + t1->mat[2][3]
	  - t2->mat[2][2]*d2 - t2->mat[2][3];

	if( AlgMatrixLUSolve(a, 3, b, 1) ){
	  rtnVtx.vtX = -1.0;
	  rtnVtx.vtY = -1.0;
	  errNum = WLZ_ERR_ALG;
	} else {
	  rtnVtx.vtX = b[0];
	  rtnVtx.vtY = b[1];
	}
	WlzFreeAffineTransform(t2);
      }
      WlzFreeAffineTransform(t1);
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnVtx;
}

/*!
* \return				The angle of intersection in radians.
* \ingroup      WlzSectionTransform
* \brief	Finds the angle of the ine of intersection of two views.
* \param	v1			First view.
* \param	v2			Second view.
* \param	dstErr			Destination pointer for an error
*					code, may be NULL.
*/
double Wlz3DViewGetIntersectionAngle(
  WlzThreeDViewStruct *v1,
  WlzThreeDViewStruct *v2,
  WlzErrorNum	*dstErr)
{
  double	rtnAngle;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzAffineTransform	*t1, *t2;
  double	vector_prod[3], l[3];
  int		i, j;

  /* check view structs - non-NULL and initialised */
  if( (v1 == NULL) || (v2 == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if( !v1->initialised || !v2->initialised ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else {
    if( t1 = WlzAffineTransformInverse(v1->trans, &errNum) ){
      if( t2 = WlzAffineTransformInverse(v2->trans, &errNum) ){
	for(i=3; i < 6; i++){
	  vector_prod[i%3] =
	    (t1->mat[(i+1)%3][2] * t2->mat[(i-1)%3][2] - 
	     t1->mat[(i-1)%3][2] * t2->mat[(i+1)%3][2]);
	}

	for(i=0; i < 3; i++ ){
	  l[i] = 0.0;
	  for(j=0; j < 3; j++){
	    l[i] += v2->trans->mat[i][j] * vector_prod[j];
	  }
	}

	/* check for undefined angle - from planes with the same normal */
	if( l[1] == 0.0 && l[0]== 0.0 )
	  rtnAngle = 0.0;
	else
	  rtnAngle = viewStructAtan2(l[1], l[0]);

	WlzFreeAffineTransform(t2);
      }
      WlzFreeAffineTransform(t1);
    }
  }
 
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnAngle;
}

/*!
* \return				The number of vertices.
* \ingroup      WlzSectionTransform
* \brief	Gets the vertices of intersection between a section and the
bounding box of the reference object. The vertices
are returned in order to be used for display etc..
* \param	viewStr			The section view.
* \param	rtnVtxs			Array of 12 vertices to return values.
* \param	dstErr			Destination pointer for an error code,
* 					may be NULL.
*/
int Wlz3DViewGetBoundingBoxIntersection(
  WlzThreeDViewStruct	*viewStr,
  WlzDVertex3		*rtnVtxs,
  WlzErrorNum		*dstErr)
{
  WlzDVertex3	bbMin, bbMax;

  /* set up bounding box */
  bbMin.vtX = viewStr->ref_obj->domain.p->kol1;
  bbMin.vtY = viewStr->ref_obj->domain.p->line1;
  bbMin.vtZ = viewStr->ref_obj->domain.p->plane1;
  bbMax.vtX = viewStr->ref_obj->domain.p->lastkl;
  bbMax.vtY = viewStr->ref_obj->domain.p->lastln;
  bbMax.vtZ = viewStr->ref_obj->domain.p->lastpl;

  return Wlz3DViewGetGivenBBIntersection(viewStr, bbMin, bbMax,
					 rtnVtxs, dstErr);
}

/*!
* \return	The number of vertices in the intersection.
* \ingroup	WlzSectionTransform
* \brief	Gets the vertices of intersection between a section and
*               the given bounding box. The vertices are returned in 
*               order to be used for display etc.
* \param	viewStr			The section View.
* \param	dstSizeArrayVtxs	Number of intersection verticies
* 					allocated, always 12 unless an
*					error occurs.
* \param	dstArrayVtxs		Destination pointer for the array for
* 					computed vertices. The array will have
* 					12 verticies allocated by AlcMalloc().
* \param	dstErr			Destination pointer for error code, may
* 					be NULL.
*/
int Wlz3DViewGetBoundingBoxIntersectionA(
  WlzThreeDViewStruct	*viewStr,
  int 			*dstSizeArrayVtxs,
  WlzDVertex3 		**dstArrayVtxs,
  WlzErrorNum 		*dstErr)
{
  int		arraySz = 0,
  		numVtxs = 0;
  WlzDVertex3	*vtxs = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(viewStr == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(viewStr->ref_obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    if( vtxs = (WlzDVertex3 *) AlcMalloc(sizeof(WlzDVertex3) * 12) ){
      numVtxs = Wlz3DViewGetBoundingBoxIntersection(viewStr,
	  vtxs, dstErr);
      arraySz = 12;
    }
    else {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(dstSizeArrayVtxs)
  {
    *dstSizeArrayVtxs = arraySz;
  }
  if(dstArrayVtxs)
  {
    *dstArrayVtxs = vtxs;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return numVtxs;
}

/* function:     Wlz3DViewGetGivenBBIntersection    */
/*! 
* \ingroup      WlzSectionTransform
* \brief        get the vertices of intersection between a section and
*		the given bounding box. The vertices
*		are returned in order to be used for display etc.
*
* \return       The number of vertices in the intersection
* \param    viewStr	The section View
* \param    bbMin	Vertex defining the minimum values of the bounding box
* \param    bbMax	Vertex defining the maximum values of the bounding box
* \param    rtnVtxs	Array of returned vertices, must be an array of at least 12
* \param    dstErr	error return
* \par      Source:
*                Wlz3DViewStructUtils.c
*/
int Wlz3DViewGetGivenBBIntersection(
  WlzThreeDViewStruct	*viewStr,
  WlzDVertex3		bbMin,
  WlzDVertex3		bbMax,
  WlzDVertex3		*rtnVtxs,
  WlzErrorNum		*dstErr)
{
  int		numVtxs=0, intersect[12], i, j=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  double	a[4], val1, val2;
  double	x1, x2, y1, y2, z1, z2;
  WlzDVertex3	vtxs[12][2], tmpVtxs[12];
  int		faces[12][2];

  /* check the view structure */
  if( viewStr == NULL ){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if( !viewStr->initialised ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  /* check the vertices - the user must supply space for at least 12 */
  if( (errNum == WLZ_ERR_NONE) && (rtnVtxs == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }

  /* calculate the plane parameters */
  Wlz3DViewGetPlaneEqn(viewStr, a + 0, a + 1, a + 2, a + 3);

  /* set up an array of lines (2 vertices each) */
  x1 = bbMin.vtX;
  x2 = bbMax.vtX;
  y1 = bbMin.vtY;
  y2 = bbMax.vtY;
  z1 = bbMin.vtZ;
  z2 = bbMax.vtZ;

  /* 3 lines from the minimum vertex values */
  vtxs[0][0].vtX = x1; vtxs[0][0].vtY = y1; vtxs[0][0].vtZ = z1;
  vtxs[0][1].vtX = x2; vtxs[0][1].vtY = y1; vtxs[0][1].vtZ = z1;
  vtxs[1][0].vtX = x1; vtxs[1][0].vtY = y1; vtxs[1][0].vtZ = z1;
  vtxs[1][1].vtX = x1; vtxs[1][1].vtY = y2; vtxs[1][1].vtZ = z1;
  vtxs[2][0].vtX = x1; vtxs[2][0].vtY = y1; vtxs[2][0].vtZ = z1;
  vtxs[2][1].vtX = x1; vtxs[2][1].vtY = y1; vtxs[2][1].vtZ = z2;
  faces[0][0] = 1; faces[0][1] = 2;
  faces[1][0] = 1; faces[1][1] = 5;
  faces[2][0] = 2; faces[2][1] = 5;

  /* 3 lines from vertex opposite to first, same plane */
  vtxs[3][0].vtX = x2; vtxs[3][0].vtY = y2; vtxs[3][0].vtZ = z1;
  vtxs[3][1].vtX = x1; vtxs[3][1].vtY = y2; vtxs[3][1].vtZ = z1;
  vtxs[4][0].vtX = x2; vtxs[4][0].vtY = y2; vtxs[4][0].vtZ = z1;
  vtxs[4][1].vtX = x2; vtxs[4][1].vtY = y1; vtxs[4][1].vtZ = z1;
  vtxs[5][0].vtX = x2; vtxs[5][0].vtY = y2; vtxs[5][0].vtZ = z1;
  vtxs[5][1].vtX = x2; vtxs[5][1].vtY = y2; vtxs[5][1].vtZ = z2;
  faces[3][0] = 1; faces[3][1] = 4;
  faces[4][0] = 1; faces[4][1] = 3;
  faces[5][0] = 3; faces[5][1] = 4;

  /* 3 lines from vertex opposite to first, same line */
  vtxs[6][0].vtX = x2; vtxs[6][0].vtY = y1; vtxs[6][0].vtZ = z2;
  vtxs[6][1].vtX = x1; vtxs[6][1].vtY = y1; vtxs[6][1].vtZ = z2;
  vtxs[7][0].vtX = x2; vtxs[7][0].vtY = y1; vtxs[7][0].vtZ = z2;
  vtxs[7][1].vtX = x2; vtxs[7][1].vtY = y2; vtxs[7][1].vtZ = z2;
  vtxs[8][0].vtX = x2; vtxs[8][0].vtY = y1; vtxs[8][0].vtZ = z2;
  vtxs[8][1].vtX = x2; vtxs[8][1].vtY = y1; vtxs[8][1].vtZ = z1;
  faces[6][0] = 2; faces[6][1] = 6;
  faces[7][0] = 3; faces[7][1] = 6;
  faces[8][0] = 2; faces[8][1] = 3;

  /* 3 lines from vertex opposite to first, same kol */
  vtxs[9][0].vtX = x1; vtxs[9][0].vtY = y2; vtxs[9][0].vtZ = z2;
  vtxs[9][1].vtX = x2; vtxs[9][1].vtY = y2; vtxs[9][1].vtZ = z2;
  vtxs[10][0].vtX = x1; vtxs[10][0].vtY = y2; vtxs[10][0].vtZ = z2;
  vtxs[10][1].vtX = x1; vtxs[10][1].vtY = y1; vtxs[10][1].vtZ = z2;
  vtxs[11][0].vtX = x1; vtxs[11][0].vtY = y2; vtxs[11][0].vtZ = z2;
  vtxs[11][1].vtX = x1; vtxs[11][1].vtY = y2; vtxs[11][1].vtZ = z1;
  faces[9][0] = 4; faces[9][1] = 6;
  faces[10][0] = 5; faces[10][1] = 6;
  faces[11][0] = 4; faces[11][1] = 5;

  /* find which lines intersect the section plane and where */
  for(i=0; i < 12; i++){
    val1 = a[0]*vtxs[i][0].vtX + a[1]*vtxs[i][0].vtY +
      a[2]*vtxs[i][0].vtZ + a[3];
    val2 = a[0]*vtxs[i][1].vtX + a[1]*vtxs[i][1].vtY +
      a[2]*vtxs[i][1].vtZ + a[3];
    if( val1 > DBL_EPSILON ){
      if( val2 > DBL_EPSILON ){
	intersect[i] = 0;
      }
      else if( val2 < DBL_EPSILON ){
	intersect[i] = 1;
	tmpVtxs[i].vtX = (-val2*vtxs[i][0].vtX + val1*vtxs[i][1].vtX)/
	  (val1 - val2);
	tmpVtxs[i].vtY = (-val2*vtxs[i][0].vtY + val1*vtxs[i][1].vtY)/
	  (val1 - val2);
	tmpVtxs[i].vtZ = (-val2*vtxs[i][0].vtZ + val1*vtxs[i][1].vtZ)/
	  (val1 - val2);
      }
      else {
	intersect[i] = 1;
	tmpVtxs[i] = vtxs[i][1];
      }
    }
    else if( val1 < -DBL_EPSILON ){
      if( val2 < -DBL_EPSILON ){
	intersect[i] = 0;
      }
      else if( val2 > DBL_EPSILON ){
	intersect[i] = 1;
	tmpVtxs[i].vtX = (val2*vtxs[i][0].vtX - val1*vtxs[i][1].vtX)/
	  (-val1 + val2);
	tmpVtxs[i].vtY = (val2*vtxs[i][0].vtY - val1*vtxs[i][1].vtY)/
	  (-val1 + val2);
	tmpVtxs[i].vtZ = (val2*vtxs[i][0].vtZ - val1*vtxs[i][1].vtZ)/
	  (-val1 + val2);
      }
      else {
	intersect[i] = 1;
	tmpVtxs[i] = vtxs[i][1];
      }
    }
    else {
      intersect[i] = 1;
      tmpVtxs[i] = vtxs[i][0];
    }
  } 
	
  /* now find correct ordering of the vertices
     - vertices are consecutive if they lie on edges of the
     same face - except if ALL vertices are on the same face - beware */
  numVtxs = 0;
  for(i=0; i < 12; i++){
    if( intersect[i] ){
      if( numVtxs ){
	/* check if consecutive to previous */
	if((faces[i][0] == faces[j][0]) || (faces[i][0] == faces[j][1]) ||
	   (faces[i][1] == faces[j][0]) || (faces[i][1] == faces[j][1])){
	  rtnVtxs[numVtxs++] = tmpVtxs[i];
	  intersect[i]=0;
	  j = i;
	  i = 0;
	}
      }
      else {
	/* first vertex */
	rtnVtxs[numVtxs++] = tmpVtxs[i];
	intersect[i]=0;
	j = i;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return numVtxs;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the fixed point coordinates from a 3D view.
* \param	vs			Given view.
* \param	dstX			Destination pointer for the fixed
*					point X coordinate.
* \param	dstY			Destination pointer for the fixed
*					point Y coordinate.
* \param	dstZ			Destination pointer for the fixed
*					point Z coordinate.
*/
WlzErrorNum Wlz3DViewGetFixed(
  WlzThreeDViewStruct	*vs,
  double		*dstX,
  double		*dstY,
  double		*dstZ)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstX = vs->fixed.vtX;
    *dstY = vs->fixed.vtY;
    *dstZ = vs->fixed.vtZ;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the fixed point coordinates in a 3D view.
* \param	vs			Given view.
* \param	x			Fixed point X coordinate.
* \param	y			Fixed point Y coordinate.
* \param	z			Fixed point Z coordinate.
*/
WlzErrorNum Wlz3DViewSetFixed(
  WlzThreeDViewStruct	*vs,
  double		x,
  double	        y,
  double		z)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->fixed.vtX = x;
    vs->fixed.vtY = y;
    vs->fixed.vtZ = z;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}


/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the angle theta from the 3D view.
* \param	vs			Given view.
* \param	dstVal			Destination pointer for theta.
*/
WlzErrorNum Wlz3DViewGetTheta(
  WlzThreeDViewStruct	*vs,
  double		*dstVal)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstVal = vs->theta;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the angle theta from the 3D view.
* \param	vs			Given view.
* \param	val			Value of theta.
*/
WlzErrorNum Wlz3DViewSetTheta(
  WlzThreeDViewStruct	*vs,
  double		val)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->theta = val;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the angle phi from the 3D view.
* \param	vs			Given view.
* \param	dstVal			Destination pointer for phi.
*/
WlzErrorNum Wlz3DViewGetPhi(
  WlzThreeDViewStruct	*vs,
  double		*dstVal)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstVal = vs->phi;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the angle phi from the 3D view.
* \param	vs			Given view.
* \param	val			Value of phi.
*/
WlzErrorNum Wlz3DViewSetPhi(
  WlzThreeDViewStruct	*vs,
  double		val)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->phi = val;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the angle zeta from the 3D view.
* \param	vs			Given view.
* \param	dstVal			Destination pointer for zeta.
*/
WlzErrorNum Wlz3DViewGetZeta(
  WlzThreeDViewStruct	*vs,
  double		*dstVal)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstVal = vs->zeta;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the angle zeta from the 3D view.
* \param	vs			Given view.
* \param	val			Value of zeta.
*/
WlzErrorNum Wlz3DViewSetZeta(
  WlzThreeDViewStruct	*vs,
  double		val)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->zeta = val;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the increment distance from the 3D view.
* \param	vs			Given view.
* \param	dstVal			Destination pointer for the
*					increment distance.
*/
WlzErrorNum Wlz3DViewGetDist(
  WlzThreeDViewStruct	*vs,
  double		*dstVal)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstVal = vs->dist;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the increment distance in the 3D view.
* \param	vs			Given view.
* \param	val			Value of the increment distance.
*/
WlzErrorNum Wlz3DViewSetDist(
  WlzThreeDViewStruct	*vs,
  double		val)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->dist = val;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the scale from the 3D view.
* \param	vs			Given view.
* \param	dstVal			Destination pointer for the
*					scale.
*/
WlzErrorNum Wlz3DViewGetScale(
  WlzThreeDViewStruct	*vs,
  double		*dstVal)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstVal = vs->scale;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the scale in the 3D view.
* \param	vs			Given view.
* \param	val			Value of the scale.
*/
WlzErrorNum Wlz3DViewSetScale(
  WlzThreeDViewStruct	*vs,
  double		val)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->scale = val;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the view mode from the 3D view.
* \param	vs			Given view.
* \param	dstVal			Destination pointer for the
*					view mode.
*/
WlzErrorNum Wlz3DViewGetViewMode(
  WlzThreeDViewStruct	*vs,
  WlzThreeDViewMode	*dstVal)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstVal = vs->view_mode;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the view mode in the 3D view.
* \param	vs			Given view.
* \param	val			Value of the view mode.
*/
WlzErrorNum Wlz3DViewSetViewMode(
  WlzThreeDViewStruct	*vs,
  WlzThreeDViewMode	val)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    switch( val ){
    case WLZ_STATUE_MODE:
    case WLZ_UP_IS_UP_MODE:
    case WLZ_FIXED_LINE_MODE:
    case WLZ_ZERO_ZETA_MODE:
    case WLZ_ZETA_MODE:
      vs->view_mode = val;
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the up vector from the 3D view.
* \param	vs			Given view.
* \param	dstX			Destination pointer for the
*					X coordinate of the up vector.
* \param	dstY			Destination pointer for the
*					Y coordinate of the up vector.
* \param	dstZ			Destination pointer for the
*					Z coordinate of the up vector.
*/
WlzErrorNum Wlz3DViewGetUp(
  WlzThreeDViewStruct	*vs,
  double		*dstX,
  double		*dstY,
  double		*dstZ)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstX = vs->up.vtX;
    *dstY = vs->up.vtY;
    *dstZ = vs->up.vtZ;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the up vector in the 3D view.
* \param	vs			Given view.
* \param	x			X coordinate of the up vector.
* \param	y			Y coordinate of the up vector.
* \param	z			Z coordinate of the up vector.
*/
WlzErrorNum Wlz3DViewSetUp(
  WlzThreeDViewStruct	*vs,
  double		x,
  double	        y,
  double		z)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->up.vtX = x;
    vs->up.vtY = y;
    vs->up.vtZ = z;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the coordinates of the second fixed point from the
*		3D view.
* \param	vs			Given view.
* \param	dstX			Destination pointer for the
*					X coordinate of the second fixed
*					point.
* \param	dstY			Destination pointer for the
*					Y coordinate of the second fixed
*					point.
* \param	dstZ			Destination pointer for the
*					Z coordinate of the second fixed
*					point.
*/
WlzErrorNum Wlz3DViewGetFixed2(
  WlzThreeDViewStruct	*vs,
  double		*dstX,
  double		*dstY,
  double		*dstZ)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstX = vs->fixed_2.vtX;
    *dstY = vs->fixed_2.vtY;
    *dstZ = vs->fixed_2.vtZ;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the coordinates of the second fixed point in the
*		3D view.
* \param	vs			Given view.
* \param	x			X coordinate of the second fixed
*					point.
* \param	y			Y coordinate of the second fixed
*					point.
* \param	z			Z coordinate of the second fixed
*					point.
*/
WlzErrorNum Wlz3DViewSetFixed2(
  WlzThreeDViewStruct	*vs,
  double		x,
  double	        y,
  double		z)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->fixed_2.vtX = x;
    vs->fixed_2.vtY = y;
    vs->fixed_2.vtZ = z;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the fixed line angle from the 3D view.
* \param	vs			Given view.
* \param	dstVal			Destination pointer for the fixed
*					line angle.
*/
WlzErrorNum Wlz3DViewGetFixedLineAngle(
  WlzThreeDViewStruct	*vs,
  double		*dstVal)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    *dstVal = vs->fixed_line_angle;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Sets the fixed line angle in the 3D view.
* \param	vs			Given view.
* \param	val			The fixed line angle.
*/
WlzErrorNum Wlz3DViewSetFixedLineAngle(
  WlzThreeDViewStruct	*vs,
  double		val)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( vs ){
    vs->fixed_line_angle = val;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the section maximum values from the 3D view.
* \param	vs			Given view.
* \param	dstX			Destination pointer for the
*					X section maximum coordinate.
* \param	dstY			Destination pointer for the
*					Y section maximum coordinate.
* \param	dstZ			Destination pointer for the
*					Z section maximum coordinate.
*/
WlzErrorNum Wlz3DViewGetMaxvals(
  WlzThreeDViewStruct	*vs,
  double		*dstX,
  double		*dstY,
  double		*dstZ)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

 /* printf("in Wlz3DViewGetMaxvals\n"); */
	 
  if( vs ){
    *dstX = vs->maxvals.vtX;
    *dstY = vs->maxvals.vtY;
    *dstZ = vs->maxvals.vtZ;
    /* printf("(x,y,z) = (%f,%f,%f)\n", *dstX, *dstY, *dstZ); */
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
    /* printf("Wlz3DViewGetMaxvals: NULL viewStruct\n"); */
  }

  return errNum;
}

/*!
* \return				Woolz error code.
* \ingroup      WlzSectionTransform
* \brief	Gets the section minimum values from the 3D view.
* \param	vs			Given view.
* \param	dstX			Destination pointer for the
*					X section maximum coordinate.
* \param	dstY			Destination pointer for the
*					Y section maximum coordinate.
* \param	dstZ			Destination pointer for the
*					Z section maximum coordinate.
*/
WlzErrorNum Wlz3DViewGetMinvals(
  WlzThreeDViewStruct	*vs,
  double		*dstX,
  double		*dstY,
  double		*dstZ)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* printf("in Wlz3DViewGetMinvals\n"); */
 	
   if( vs ){
    *dstX = vs->minvals.vtX;
    *dstY = vs->minvals.vtY;
    *dstZ = vs->minvals.vtZ;
  }
  else {
    errNum = WLZ_ERR_PARAM_NULL;
  }

  return errNum;
}

/*!
* \return	void
* \brief	Computes the parameters of the equation of the plane
*		defined by the given 3D view.
*
*		Computes the parameters of the equation of the plane
*		\f[ax + by + cz + d = 0\f] defined by the given 3D view.
*		None of the destination pointers may be NULL and the
*		3D view is assumed valid although only its theta, phi
*		distance increment and fixed point are used.
* \param	view			Given view.
* \param	dstA			Destination pointer for the
*					X parameter.
* \param	dstB			Destination pointer for the
*					Y parameter.
* \param	dstC			Destination pointer for the
*					Z parameter.
* \param	dstD			Destination pointer for the
*					other parameter.
*/
void		Wlz3DViewGetPlaneEqn(WlzThreeDViewStruct *view,
				     double *dstA, double *dstB,
				     double *dstC, double *dstD)
{
  double	a,
  		b,
		c,
		d,
		sPhi;

  sPhi = sin(view->phi);
  a = sPhi * cos(view->theta);
  b = sPhi * sin(view->theta);
  c = cos(view->phi);
  d = -(a * (view->fixed.vtX + (a * view->dist)) +
	b * (view->fixed.vtY + (b * view->dist)) +
	c * (view->fixed.vtZ + (c * view->dist)));
  *dstA = a;
  *dstB = b;
  *dstC = c;
  *dstD = d;
}

/*!
* \return				Non zero if there is an intersection
*					otherwise zero.
* \brief	Tests for an intersection between the plane defined by the
* 		given 3D view and the given axis aligned bounding box (AABB).
* \param	view			Given view which defines a plane.
* \param	box			Given axis aligned bounding box.
*/
int		Wlz3DViewIntersectAABB(WlzThreeDViewStruct *view,
				       WlzDBox3 box)
{
  int		intersect = 0;
  double	p[4];

  /* Compute the eqn of the plane. */
  Wlz3DViewGetPlaneEqn(view, p + 0, p + 1, p + 2, p + 3);
  /* Check for intersection. */
  intersect = WlzGeomPlaneAABBIntersect(p[0], p[1], p[2], p[3], box);
  return(intersect);
}
