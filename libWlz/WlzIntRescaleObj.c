#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzIntRescaleObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzIntRescaleObj.c
* \author       Richard Baldock, Bill Hill
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
* \brief	Rescales a Woolz object using an integral scale.
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

static int 			check_intvs(
				  WlzInterval *intvs,
				  int nintvs);
static WlzObject 		*WlzIntRescaleObj2D(
				  WlzObject *obj,
				  int scale,
				  int expand,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzIntRescaleObj3D(
				  WlzObject *obj,
				  int scale,
				  int expand,
				  WlzErrorNum *dstErr);

/*!
* \return	Rescaled object.
* \ingroup	WlzTransform
* \brief	Rescales the given object using an integer scale.
* \param	gObj			Given object.
* \param	scale			Integer scale factor.
* \param	expand			If zero use \f$\frac{1}{scale}\f$.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzIntRescaleObj(WlzObject *gObj, int scale, int expand,
				WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(scale < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	if(gObj->domain.i == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(gObj->domain.i->type == WLZ_EMPTY_DOMAIN)
	{
	  rObj = WlzMakeEmpty(&errNum);
	}
	else if(scale == 1)
	{
          rObj = WlzMakeMain(gObj->type, gObj->domain, gObj->values,
			     NULL, NULL, &errNum);
	}
	else
	{
	  rObj = (gObj->type == WLZ_2D_DOMAINOBJ)?
	         WlzIntRescaleObj2D(gObj, scale, expand, &errNum):
	         WlzIntRescaleObj3D(gObj, scale, expand, &errNum);
	}
	break;
      case WLZ_TRANS_OBJ: /* FALLTHROUGH */
      case WLZ_2D_POLYGON: /* FALLTHROUGH */
      case WLZ_BOUNDLIST: /* FALLTHROUGH */
      case WLZ_3D_POLYGON:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      case WLZ_EMPTY_OBJ:
	rObj = WlzMakeEmpty(&errNum);
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Number of intervals.
* \ingroup	WlzTransform
* \brief	Check and count intervals.
* \param	intvs		Given intervals.
* \param	nintvs		Number of intervals.
*/
static int check_intvs(
  WlzInterval	*intvs,
  int		nintvs)
{
  int	i;

  /* stopping condition */
  if( nintvs == 0 )
    return( 0 );

  /* check tail of the list */
  nintvs = check_intvs(intvs+1, nintvs-1) + 1;

  /* check this interval */
  if( intvs->iright < intvs->ileft ){
    intvs->iright = intvs->ileft;
  }

  /* return if one interval */
  if( nintvs == 1 )
    return( nintvs );

  /* check next interval */
  if( (intvs[1].ileft - intvs[0].iright) < 2 ){ /* merge condition */
    intvs[0].iright = intvs[1].iright;
    for(i=2; i < nintvs; i++)
      intvs[i-1] = intvs[i];
    return( nintvs - 1 );
  }

  return( nintvs );
}

/*!
* \return	Rescaled object.
* \ingroup	WlzTransform
* \brief	Rescales the given 2D domain object using an integer scale.
* \param	obj			Given object.
* \param	scale			Integer scale factor.
* \param	expand			If zero use \f$\frac{1}{scale}\f$.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzIntRescaleObj2D(
  WlzObject	*obj,
  int		scale,
  int		expand,
  WlzErrorNum	*dstErr)
{
  WlzObject		*rtnObj=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzInterval		*intvls;
  int			k1, kl, l1, ll, l, num_intvls;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check expand or contract */
  if( expand )
  {
    k1 = obj->domain.i->kol1   * scale;
    kl = obj->domain.i->lastkl * scale + scale - 1;
    l1 = obj->domain.i->line1  * scale;
    ll = obj->domain.i->lastln * scale + scale - 1;
  }
  else {
    k1 = obj->domain.i->kol1   / scale;
    kl = obj->domain.i->lastkl / scale;
    l1 = obj->domain.i->line1  / scale;
    ll = obj->domain.i->lastln / scale;
  }

  /* create a new object */
  if( domain.i = WlzMakeIntervalDomain(obj->domain.i->type,
				       l1, ll, k1, kl, &errNum) ){
    values.core = NULL;
    rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, 
			 NULL, NULL, NULL);
  }

  /* fill in the intervals */
  if( errNum == WLZ_ERR_NONE){
    if( domain.i->type == WLZ_INTERVALDOMAIN_INTVL )
    {
      int intvline_offset, max_offset;
      WlzIntervalLine *intvline;

      num_intvls = WlzIntervalCount(obj->domain.i, NULL);
      num_intvls = expand ? num_intvls * scale : num_intvls;
      intvls = (WlzInterval *)AlcMalloc(sizeof(WlzInterval) * num_intvls);
      domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)intvls,
      					   NULL);
      max_offset = obj->domain.i->lastln - obj->domain.i->line1;

      for(l=l1; l <= ll; l++)
      {
	int		i;

	intvline_offset = (expand?l/scale:l*scale) - obj->domain.i->line1;
	intvline_offset = WLZ_MAX(intvline_offset, 0);
	intvline_offset = WLZ_MIN(intvline_offset, max_offset);
	intvline = obj->domain.i->intvlines + intvline_offset;

	for(i=0; i < intvline->nintvs; i++)
	{
	  intvls[i].ileft  = (intvline->intvs + i)->ileft;
	  intvls[i].iright = (intvline->intvs + i)->iright;
	  if( expand )
	  {
	    intvls[i].ileft  *= scale;
	    intvls[i].iright *= scale;
	    intvls[i].iright += scale - 1;
	  }
	  else
	  {
	    intvls[i].ileft  /= scale;
	    intvls[i].iright /= scale;
	  }
	}

	i = check_intvs(intvls, i);
	WlzMakeInterval(l, domain.i, i, intvls);
	intvls += i;
      }
      (void) WlzStandardIntervalDomain( domain.i );
    }
  }

  /* create the valuetable */
  if( (errNum == WLZ_ERR_NONE) && obj->values.core )
  {
    WlzIntervalWSpace	iwsp;
    WlzGreyWSpace	gwsp;
    WlzPixelV 		backgrnd;
    WlzGreyValueWSpace	*gVWSp = NULL;
    WlzGreyType		gtype;

    backgrnd = WlzGetBackground(obj, NULL);
    if( values.v = WlzNewValueTb(rtnObj, obj->values.v->type,
				 backgrnd, &errNum) ){

      rtnObj->values = WlzAssignValues(values, NULL);
      
      /* fill in the grey-values */
      WlzInitGreyScan(rtnObj, &iwsp, &gwsp);
      gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
      gtype = WlzGreyTableTypeToGreyType(obj->values.v->type, NULL);
      while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE )
      {
	int k;
	int lp = expand ? iwsp.linpos/scale : iwsp.linpos*scale;

	for( k=0; k <= (iwsp.rgtpos - iwsp.lftpos); k++ )
	{
	  int	kp = expand ? (k+iwsp.lftpos)/scale : (k+iwsp.lftpos)*scale;

	  WlzGreyValueGet(gVWSp, 0, (double) lp, (double) kp);

	  switch(gtype)
	  {
	  case WLZ_GREY_INT: 
	    gwsp.u_grintptr.inp[k] = (*(gVWSp->gVal)).inv;
	    break;
	  case WLZ_GREY_SHORT:
	    gwsp.u_grintptr.shp[k] = (*(gVWSp->gVal)).shv;
	    break;
	  case WLZ_GREY_UBYTE:
	    gwsp.u_grintptr.ubp[k] = (*(gVWSp->gVal)).ubv;
	    break;
	  case WLZ_GREY_FLOAT:
	    gwsp.u_grintptr.flp[k] = (*(gVWSp->gVal)).flv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    gwsp.u_grintptr.dbp[k] = (*(gVWSp->gVal)).dbv;
	    break;
	  case WLZ_GREY_RGBA:
	    gwsp.u_grintptr.rgbp[k] = (*(gVWSp->gVal)).rgbv;
	    break;
	  }
	}
      }
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
      WlzGreyValueFreeWSp(gVWSp);
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

/*!
* \return	Rescaled object.
* \ingroup	WlzTransform
* \brief	Rescales the given 3D domain object using an integer scale.
* \param	gObj			Given object.
* \param	scale			Integer scale factor.
* \param	expand			If zero use \f$\frac{1}{scale}\f$.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzIntRescaleObj3D(WlzObject *gObj, int scale, int expand,
				    WlzErrorNum *dstErr)
{
  int		gPIdx,
  		nPIdx;
  WlzDomain	gDom,
  		nDom;
  WlzValues	gVal,
  		nVal,
		dumVal;
  WlzObject	*gTObj = NULL,
  		*nTObj = NULL,
		*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIBox3	nBox;

  nDom.core = NULL;
  nVal.core = NULL;
  dumVal.core = NULL;
  gDom = gObj->domain;
  gVal = gObj->values;
  if(expand)
  {
    nBox.xMin = gDom.p->kol1 * scale;
    nBox.xMax = ((gDom.p->lastkl + 1) * scale) - 1;
    nBox.yMin = gDom.p->line1 * scale;
    nBox.yMax = ((gDom.p->lastln + 1) * scale) - 1;
    nBox.zMin = gDom.p->plane1 * scale;
    nBox.zMax = ((gDom.p->lastpl + 1) * scale) - 1;
  }
  else
  {
    nBox.xMin = gDom.p->kol1 / scale;
    nBox.xMax = gDom.p->lastkl / scale;
    nBox.yMin = gDom.p->line1 / scale;
    nBox.yMax = gDom.p->lastln / scale;
    nBox.zMin = gDom.p->plane1 / scale;
    nBox.zMax = gDom.p->lastpl / scale;
  }
  nDom.p = WlzMakePlaneDomain(gDom.p->type, nBox.zMin, nBox.zMax,
			      nBox.yMin, nBox.yMax, nBox.xMin, nBox.xMax,
			      &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    if(gVal.core && (gVal.core->type != WLZ_EMPTY_OBJ))
    {
      nVal.vox = WlzMakeVoxelValueTb(gVal.vox->type, nBox.zMin, nBox.zMax,
				     WlzGetBackground(gObj, NULL),
				     NULL, &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nPIdx = nBox.zMin;
    while((errNum == WLZ_ERR_NONE) && (nPIdx <= nBox.zMax))
    {
      gPIdx = (expand)? nPIdx / scale: nPIdx * scale;
      if(nVal.vox)
      {
        gTObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
                            *(gDom.p->domains + gPIdx),
                            *(gVal.vox->values + gPIdx),
                            NULL, NULL, &errNum);
      }
      else
      {
        gTObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
                            *(gDom.p->domains + gPIdx),
                            dumVal,
                            NULL, NULL, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	nTObj = WlzIntRescaleObj2D(gTObj, scale, expand, &errNum);
      }
      (void )WlzFreeObj(gTObj);
      gTObj = NULL;
      if(errNum == WLZ_ERR_NONE)
      {
        *(nDom.p->domains + nPIdx) = WlzAssignDomain(nTObj->domain, NULL);
        if(nVal.vox)
        {
          *(nVal.vox->values + nPIdx) = WlzAssignValues(nTObj->values, NULL);
        }
      }
      (void )WlzFreeObj(nTObj);
      nTObj = NULL;
      ++nPIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(gObj->type, nDom, nVal, NULL, NULL, &errNum);
  }
  else
  {
    (void )WlzFreePlaneDomain(nDom.p);
    (void )WlzFreeVoxelValueTb(nVal.vox);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}
