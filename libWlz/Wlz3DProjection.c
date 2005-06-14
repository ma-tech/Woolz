#pragma ident "MRC HGU $Id$"
/*!
* \file         Wlz3DProjection.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Mon Mar  7 14:44:45 2005
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzTransform
* \brief        Generate the projection or back-projection of a
*		woolz domain object from 3D to 2D and visa versa.
*               
* \todo         Implement grey-value integration function call.
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/* function:     WlzGetProjectionFromObject    */
/*! 
* \ingroup      WlzTransform
* \brief        Use the view transform to define a projection from
*		3D to 2D. Currently only the domain is projected as
*		an opaque shadow.
*
* \return       projection object
* \param    obj	source 3D object
* \param    viewStr	view structure defining the projection
* \param    intFunc	grey-value summation function
* \param    intFuncData data to be passed to the integration function
* \param    dstErr	error return
* \par      Source:
*                Wlz3DProjection.c
*/
WlzObject *WlzGetProjectionFromObject(
  WlzObject		*obj,
  WlzThreeDViewStruct 	*viewStr,
  Wlz3DProjectionIntFn 	intFunc,
  void			*intFuncData,
  WlzErrorNum		*dstErr)
{
  WlzObject		*rtnObj=NULL,
  			*obj1;
  WlzThreeDViewStruct	*viewStr1=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzGreyType		srcGType, dstGType;
  WlzPixelV		pixval;
  WlzPixelP		pixptr;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyValueWSpace	*gVWSp = NULL;
  WlzDVertex3		vtx, vtx1;
  double		x, y, z;
  double		*s_to_x=NULL;
  double		*s_to_y=NULL;
  double		*s_to_z=NULL;
  int			k, xp, yp, p, s, sp;
  int			length, size, occupiedFlg;
  WlzErrorNum 	errNum=WLZ_ERR_NONE;

  /* check inputs */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if( obj->type != WLZ_3D_DOMAINOBJ ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  if( (errNum == WLZ_ERR_NONE) && (viewStr == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* create new view transform */
  if( errNum == WLZ_ERR_NONE ){
    if( viewStr1 = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum) ){
      /* need to worry about fixed line mode here sometime */
      viewStr1->fixed = viewStr->fixed;
      viewStr1->theta = viewStr->theta;
      viewStr1->phi = viewStr->phi;
      viewStr1->zeta = viewStr->zeta;
      viewStr1->dist = viewStr->dist;
      viewStr1->scale = viewStr->scale;
      viewStr1->voxelSize[0] = viewStr->voxelSize[0];
      viewStr1->voxelSize[1] = viewStr->voxelSize[1];
      viewStr1->voxelSize[2] = viewStr->voxelSize[2];
      viewStr1->voxelRescaleFlg = viewStr->voxelRescaleFlg;
      viewStr1->interp = viewStr->interp;
      viewStr1->view_mode = viewStr->view_mode;
      viewStr1->up = viewStr->up;

      /* now intialize it */
      /* could optimise by setting fixed point to object centre */
      if( (errNum = WlzInit3DViewStruct(viewStr1, obj)) != WLZ_ERR_NONE ){
	WlzFree3DViewStruct(viewStr1);
	viewStr1 = NULL;
      }
    }
  }

  /* set up orthogonal line parameters & luts */
  if( errNum == WLZ_ERR_NONE ){
    length = WLZ_NINT(viewStr1->maxvals.vtZ) -
      WLZ_NINT(viewStr1->minvals.vtZ) + 1;
    s_to_x = (double *) AlcMalloc(sizeof(double) * length );
    s_to_y = (double *) AlcMalloc(sizeof(double) * length );
    s_to_z = (double *) AlcMalloc(sizeof(double) * length );

    /* transform a perpendicular vector */
    vtx.vtX = 0.0;
    vtx.vtY = 0.0;
    vtx.vtZ = 1.0;
    Wlz3DSectionTransformInvVtx(&vtx, viewStr1);
    vtx1.vtX = 0.0;
    vtx1.vtY = 0.0;
    vtx1.vtZ = 0.0;
    Wlz3DSectionTransformInvVtx(&vtx1, viewStr1);
    vtx.vtX -= vtx1.vtX;
    vtx.vtY -= vtx1.vtY;
    vtx.vtZ -= vtx1.vtZ;

    /* assign lut values */
    s = WLZ_NINT(viewStr1->minvals.vtZ) - viewStr1->dist;
    for(sp=0; sp < length; sp++, s++){
      s_to_x[sp] = s * vtx.vtX;
      s_to_y[sp] = s * vtx.vtY;
      s_to_z[sp] = s * vtx.vtZ;
    }
  }

  /* if there is an integration function then allocate space for
     the grey-level array */
  if( (errNum == WLZ_ERR_NONE) && (intFunc) ){
    srcGType = WlzGreyTypeFromObj(obj, &errNum);
    switch( srcGType ){
    case WLZ_GREY_LONG:
      size = sizeof(long)*length;
      break;
    case WLZ_GREY_INT:
      size = sizeof(int)*length;
      break;
    case WLZ_GREY_SHORT:
      size = sizeof(short)*length;
      break;
    case WLZ_GREY_UBYTE:
      size = sizeof(UBYTE)*length;
      break;
    case WLZ_GREY_FLOAT:
      size = sizeof(float)*length;
      break;
    case WLZ_GREY_DOUBLE:
      size = sizeof(double)*length;
      break;
    case WLZ_GREY_RGBA:
      size = sizeof(int)*length;
      break;
    }
    if( (pixptr.p.inp = (int *) AlcMalloc(size)) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    pixptr.type = srcGType;

    /* set up the grey-value workspace for random access */
    gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
  }

  /* create rectangular projection image */
  if( errNum == WLZ_ERR_NONE ){
    if( domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				     WLZ_NINT(viewStr1->minvals.vtY),
				     WLZ_NINT(viewStr1->maxvals.vtY),
				     WLZ_NINT(viewStr1->minvals.vtX),
				     WLZ_NINT(viewStr1->maxvals.vtX),
					 &errNum) ){
      values.core = NULL;
      if( rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL,
			       &errNum) ){
	/* note the grey-values required are determined by the integration
	   function. Here we use UBYTE and reset later if needed */
	dstGType = WLZ_GREY_UBYTE;
	pixval.type = WLZ_GREY_UBYTE;
	pixval.v.ubv = (UBYTE) 0;
	if( values.v = WlzNewValueTb(rtnObj,
				     WlzGreyTableType(WLZ_GREY_TAB_RECT,
						      dstGType, NULL),
				     pixval, &errNum) ){
	  rtnObj->values = WlzAssignValues(values, &errNum);
	}
	else {
	  WlzFreeObj(rtnObj);
	  rtnObj = NULL;
	}
      }
      else {
	WlzFreeDomain(domain);
	domain.core = NULL;
      }
    }
  }

  /* scan image setting values */
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(rtnObj, &iwsp, &gwsp);
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
      yp = iwsp.linpos - WLZ_NINT(viewStr1->minvals.vtY);
      for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++){
	xp = k - WLZ_NINT(viewStr1->minvals.vtX);
	vtx.vtX = viewStr1->xp_to_x[xp] + viewStr1->yp_to_x[yp];
	vtx.vtY = viewStr1->xp_to_y[xp] + viewStr1->yp_to_y[yp];
	vtx.vtZ = viewStr1->xp_to_z[xp] + viewStr1->yp_to_z[yp];

	/* get the projection values */
	/* if no function then just check for occupancy */
	if( intFunc == NULL ){
	  occupiedFlg = 0;
	  sp = viewStr1->dist - WLZ_NINT(viewStr1->minvals.vtZ);
	  for(; !occupiedFlg && (sp < length); sp++){
	    x = vtx.vtX + s_to_x[sp];
	    y = vtx.vtY + s_to_y[sp];
	    z = vtx.vtZ + s_to_z[sp];
	    if( WlzInsideDomain(obj, z, y, x, &errNum) ){
	      occupiedFlg = 1;
	    }
	  }
	  sp = viewStr1->dist - WLZ_NINT(viewStr1->minvals.vtZ) - 1;
	  for(; !occupiedFlg && (sp >= 0); sp--){
	    x = vtx.vtX + s_to_x[sp];
	    y = vtx.vtY + s_to_y[sp];
	    z = vtx.vtZ + s_to_z[sp];
	    if( WlzInsideDomain(obj, z, y, x, &errNum) ){
	      occupiedFlg = 1;
	    }
	  }

	  /* set the integrated value - only UBYTE at the moment */
	  *(gwsp.u_grintptr.ubp) = occupiedFlg;
	  gwsp.u_grintptr.ubp++;
	}
	/* use integration function */
	else {
	  /* set array of pixel values */
	  for(sp=0; sp < length; sp++){
	    x = vtx.vtX + s_to_x[sp];
	    y = vtx.vtY + s_to_y[sp];
	    z = vtx.vtZ + s_to_z[sp];
	    WlzGreyValueGet(gVWSp, WLZ_NINT(z), WLZ_NINT(y),
			    WLZ_NINT(x));
	    switch( srcGType ){
	    case WLZ_GREY_LONG:
	      pixptr.p.lnp[sp] = gVWSp->gVal[0].lnv;
	      break;
	    case WLZ_GREY_INT:
	      pixptr.p.inp[sp] = gVWSp->gVal[0].inv;
	      break;
	    case WLZ_GREY_SHORT:
	      pixptr.p.shp[sp] = gVWSp->gVal[0].shv;
	      break;
	    case WLZ_GREY_UBYTE:
	      pixptr.p.ubp[sp] = gVWSp->gVal[0].ubv;
	      break;
	    case WLZ_GREY_FLOAT:
	      pixptr.p.flp[sp] = gVWSp->gVal[0].flv;
	      break;
	    case WLZ_GREY_DOUBLE:
	      pixptr.p.dbp[sp] = gVWSp->gVal[0].dbv;
	      break;
	    case WLZ_GREY_RGBA:
	      pixptr.p.rgbp[sp] = gVWSp->gVal[0].rgbv;
	      break;
	    }
	  }
	  /* call integration function and seet value */
	  intFunc(pixptr, length, viewStr1->dist - WLZ_NINT(viewStr1->minvals.vtZ),
		  intFuncData, &errNum);
	}
      }
    if(errNum == WLZ_ERR_EOO)	   /* Reset error from end of intervals */ 
    {
      errNum = WLZ_ERR_NONE;
    }

    /* if no integration function then threshold - binary only */
    if( intFunc == NULL ){
      pixval.v.ubv = 1;
      rtnObj = WlzAssignObject(rtnObj, NULL);
      if( obj1 = WlzThreshold(rtnObj, pixval, WLZ_THRESH_HIGH, &errNum) ){
	WlzFreeObj(rtnObj);
	values.core = NULL;
	rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj1->domain, values, NULL, NULL,
			     &errNum);
	WlzFreeObj(obj1);
      }
      else {
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
    }
    }
  }

  /* clear space */
  if( viewStr1 ){
    errNum = WlzFree3DViewStruct(viewStr1);
  }
  if( s_to_x ){
    AlcFree( s_to_x );
  }
  if( s_to_y ){
    AlcFree( s_to_y );
  }
  if( s_to_z ){
    AlcFree( s_to_z );
  }

  /* check error and return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}