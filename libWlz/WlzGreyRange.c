#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyRange_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreyRange.c
* \author       Richard Baldock
* \date         September 2005
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
* \brief	Computes the grey range of an object.
* 		Note that the colour range returned is the independent max
* 		and min values for each channel.
* \ingroup	WlzValuesFilters
*/

#include <stdlib.h>
#include <Wlz.h>

/* function:     WlzGreyRange    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        compute grey-range of a pixel/voxel object.
*
* \return       Woolz error number: WLZ_ERR_NONE, WLZ_ERR_OBJECT_NULL,
 WLZ_ERR_DOMAIN_NULL, WLZ_ERR_VALUES_NULL, WLZ_ERR_GREY_TYPE, 
WLZ_ERR_VOXELVALUES_TYPE, WLZ_ERR_PLANEDOMAIN_TYPE, WLZ_ERR_OBJECT_TYPE
* \param    obj	grey-level input object
* \param    min	minimum values return
* \param    max	maximum values return.
* \par      Source:
*                WlzGreyRange.c
*/
WlzErrorNum WlzGreyRange(WlzObject	*obj,
			 WlzPixelV	*min,
			 WlzPixelV	*max)
{
  WlzGreyV		v;
  WlzGreyP		g;
  int			i, nplanes, init_flag;
  WlzPixelV		lmin, lmax, Min, Max;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzObject		tempobj;
  WlzPlaneDomain	*planedm;
  WlzValues 		*values;
  WlzDomain		*domains;
  WlzGreyType		gType;
  WlzErrorNum		errNum;

  /* check for NULL object */
  if( obj == NULL ){
    return( WLZ_ERR_OBJECT_NULL );
  }

  /* check for NULL domain or values */
  if( obj->domain.core == NULL ){
    return( WLZ_ERR_DOMAIN_NULL );
  }
  if( obj->values.core == NULL ){
    return( WLZ_ERR_VALUES_NULL );
  }
    
  init_flag = 0;
  switch( obj->type ){

  case WLZ_2D_DOMAINOBJ:
    WlzInitGreyScan(obj, &iwsp, &gwsp);
    lmin.type = gwsp.pixeltype;
    lmax.type = gwsp.pixeltype;
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
      g = gwsp.u_grintptr;
      switch( gwsp.pixeltype ){

      case WLZ_GREY_INT:
	for(i=0; i<iwsp.colrmn; i++){
	  v.inv = *g.inp++;
	  if( !init_flag ){
	    lmin.v.inv = v.inv;
	    lmax.v.inv = v.inv;
	    init_flag = 1;
	  }
	  else if( v.inv > lmax.v.inv ){
	    lmax.v.inv = v.inv;
	  }
	  else if( v.inv < lmin.v.inv ){
	    lmin.v.inv = v.inv;
	  }
	}
	break;

      case WLZ_GREY_SHORT:
	for(i=0; i<iwsp.colrmn; i++){
	  v.shv = *g.shp++;
	  if( !init_flag ){
	    lmin.v.shv = v.shv;
	    lmax.v.shv = v.shv;
	    init_flag = 1;
	  }
	  else if( v.shv > lmax.v.shv ){
	    lmax.v.shv = v.shv;
	  }
	  else if( v.shv < lmin.v.shv ){
	    lmin.v.shv = v.shv;
	  }
	}
	break;

      case WLZ_GREY_UBYTE:
	for(i=0; i<iwsp.colrmn; i++){
	  v.ubv = *g.ubp++;
	  if( !init_flag ){
	    lmin.v.ubv = v.ubv;
	    lmax.v.ubv = v.ubv;
	    init_flag = 1;
	  }
	  else if( v.ubv > lmax.v.ubv ){
	    lmax.v.ubv = v.ubv;
	  }
	  else if( v.ubv < lmin.v.ubv ){
	    lmin.v.ubv = v.ubv;
	  }
	}
	break;

      case WLZ_GREY_FLOAT:
	for(i=0; i<iwsp.colrmn; i++){
	  v.flv = *g.flp++;
	  if( !init_flag ){
	    lmin.v.flv = v.flv;
	    lmax.v.flv = v.flv;
	    init_flag = 1;
	  }
	  else if( v.flv > lmax.v.flv ){
	    lmax.v.flv = v.flv;
	  }
	  else if( v.flv < lmin.v.flv ){
	    lmin.v.flv = v.flv;
	  }
	}
	break;

      case WLZ_GREY_DOUBLE:
	for(i=0; i<iwsp.colrmn; i++){
	  v.dbv = *g.dbp++;
	  if( !init_flag ){
	    lmin.v.dbv = v.dbv;
	    lmax.v.dbv = v.dbv;
	    init_flag = 1;
	  }
	  else if( v.dbv > lmax.v.dbv ){
	    lmax.v.dbv = v.dbv;
	  }
	  else if( v.dbv < lmin.v.dbv ){
	    lmin.v.dbv = v.dbv;
	  }
	}
	break;

      case WLZ_GREY_RGBA:
	for(i=0; i<iwsp.colrmn; i++){
	  v.rgbv = *g.rgbp++;
	  if( !init_flag ){
	    lmin.v.rgbv = v.rgbv;
	    lmax.v.rgbv = v.rgbv;
	    init_flag = 1;
	  }
	  else {
	    /* red */
	    if( (v.rgbv&0xff) > (lmax.v.rgbv&0xff) ){
	      lmax.v.rgbv = (lmax.v.rgbv&0xffffff00) + (v.rgbv&0xff);
	    }
	    else if( (v.rgbv&0xff) < (lmin.v.rgbv&0xff) ){
	      lmin.v.rgbv = (lmin.v.rgbv&0xffffff00) + (v.rgbv&0xff);
	    }

	    /* green */
	    if( (v.rgbv&0xff00) > (lmax.v.rgbv&0xff00) ){
	      lmax.v.rgbv = (lmax.v.rgbv&0xffff00ff) + (v.rgbv&0xff00);
	    }
	    else if( (v.rgbv&0xff00) < (lmin.v.rgbv&0xff00) ){
	      lmin.v.rgbv = (lmin.v.rgbv&0xffff00ff) + (v.rgbv&0xff00);
	    }

	    /* blue */
	    if( (v.rgbv&0xff0000) > (lmax.v.rgbv&0xff0000) ){
	      lmax.v.rgbv = (lmax.v.rgbv&0xff00ffff) + (v.rgbv&0xff0000);
	    }
	    else if( (v.rgbv&0xff0000) < (lmin.v.rgbv&0xff0000) ){
	      lmin.v.rgbv = (lmin.v.rgbv&0xff00ffff) + (v.rgbv&0xff0000);
	    }

	    /* alpha */
	    if( (v.rgbv&0xff000000) > (lmax.v.rgbv&0xff000000) ){
	      lmax.v.rgbv = (lmax.v.rgbv&0xffffff00) + (v.rgbv&0xff000000);
	    }
	    else if( (v.rgbv&0xff000000) < (lmin.v.rgbv&0xff000000) ){
	      lmin.v.rgbv = (lmin.v.rgbv&0x00ffffff) + (v.rgbv&0xff000000);
	    }
	  }
	}
	break;

      default:
	return( WLZ_ERR_GREY_TYPE );

      }
    }
    if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }
    *min = lmin;
    *max = lmax;
    break;

  case WLZ_3D_DOMAINOBJ:
    planedm = obj->domain.p;
    switch( planedm->type ){

    case WLZ_PLANEDOMAIN_DOMAIN:
      nplanes = planedm->lastpl - planedm->plane1 + 1;
      domains = planedm->domains;
      if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	return( WLZ_ERR_VOXELVALUES_TYPE );
      }
      values = obj->values.vox->values;
      tempobj.type = WLZ_2D_DOMAINOBJ;
      tempobj.plist = NULL;
      tempobj.assoc = NULL;
      for(i=0; i < nplanes; i++, domains++, values++){
	if( (*domains).core == NULL || (*values).core == NULL ){
	  continue;
	}

	tempobj.domain = *domains;
	tempobj.values = *values;
	errNum = WlzGreyRange(&tempobj, &Min, &Max);
	if( errNum != WLZ_ERR_NONE ){
	  return( errNum );
	}
	if( !init_flag ){
	  lmin = Min;
	  lmax = Max;
	  init_flag = 1;
	  continue;
	}
	
	gType = WlzGreyTableTypeToGreyType((*values).v->type, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  return(errNum);
	}
	switch(gType){

	case WLZ_GREY_INT:
	  if( Min.v.inv < lmin.v.inv ){
	    lmin.v.inv = Min.v.inv;
	  }
	  if( Max.v.inv > lmax.v.inv ){
	    lmax.v.inv = Max.v.inv;
	  }
	  break;

	case WLZ_GREY_SHORT:
	  if( Min.v.shv < lmin.v.shv ){
	    lmin.v.shv = Min.v.shv;
	  }
	  if( Max.v.shv > lmax.v.shv ){
	    lmax.v.shv = Max.v.shv;
	  }
	  break;

	case WLZ_GREY_UBYTE:
	  if( Min.v.ubv < lmin.v.ubv ){
	    lmin.v.ubv = Min.v.ubv;
	  }
	  if( Max.v.ubv > lmax.v.ubv ){
	    lmax.v.ubv = Max.v.ubv;
	  }
	  break;

	case WLZ_GREY_FLOAT:
	  if( Min.v.flv < lmin.v.flv ){
	    lmin.v.flv = Min.v.flv;
	  }
	  if( Max.v.flv > lmax.v.flv ){
	    lmax.v.flv = Max.v.flv;
	  }
	  break;

	case WLZ_GREY_DOUBLE:
	  if( Min.v.dbv < lmin.v.dbv ){
	    lmin.v.dbv = Min.v.dbv;
	  }
	  if( Max.v.dbv > lmax.v.dbv ){
	    lmax.v.dbv = Max.v.dbv;
	  }
	  break;

	case WLZ_GREY_RGBA:
	  /* red */
	  if( (Min.v.rgbv&0xff) < (lmin.v.rgbv&0xff) ){
	    lmin.v.rgbv &= ~0xff;
	    lmin.v.rgbv |= (Min.v.rgbv&0xff);
	  }
	  if( (Max.v.rgbv&0xff) > (lmax.v.rgbv&0xff) ){
	    lmin.v.rgbv &= ~0xff;
	    lmin.v.rgbv |= (Max.v.rgbv&0xff);
	  }

	  /* green */
	  if( (Min.v.rgbv&0xff00) < (lmin.v.rgbv&0xff00) ){
	    lmin.v.rgbv &= ~0xff00;
	    lmin.v.rgbv |= (Min.v.rgbv&0xff00);
	  }
	  if( (Max.v.rgbv&0xff00) > (lmax.v.rgbv&0xff00) ){
	    lmin.v.rgbv &= ~0xff00;
	    lmin.v.rgbv |= (Max.v.rgbv&0xff00);
	  }

	  /* blue */
	  if( (Min.v.rgbv&0xff0000) < (lmin.v.rgbv&0xff0000) ){
	    lmin.v.rgbv &= ~0xff0000;
	    lmin.v.rgbv |= (Min.v.rgbv&0xff0000);
	  }
	  if( (Max.v.rgbv&0xff0000) > (lmax.v.rgbv&0xff0000) ){
	    lmin.v.rgbv &= ~0xff0000;
	    lmin.v.rgbv |= (Max.v.rgbv&0xff0000);
	  }

	  /* alpha */
	  if( (Min.v.rgbv&0xff000000) < (lmin.v.rgbv&0xff000000) ){
	    lmin.v.rgbv &= ~0xff000000;
	    lmin.v.rgbv |= (Min.v.rgbv&0xff000000);
	  }
	  if( (Max.v.rgbv&0xff000000) > (lmax.v.rgbv&0xff000000) ){
	    lmin.v.rgbv &= ~0xff000000;
	    lmin.v.rgbv |= (Max.v.rgbv&0xff000000);
	  }

	  break;

	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  return(errNum);
	  /* break; */
	}
      }
      break;

    default:
      return( WLZ_ERR_PLANEDOMAIN_TYPE );

    }
    *min = lmin;
    *max = lmax;
    break;

  case WLZ_TRANS_OBJ:
    return( WlzGreyRange(obj->values.obj, min, max) );

  case WLZ_EMPTY_OBJ:
    return WLZ_ERR_NONE;

  default:
    return( WLZ_ERR_OBJECT_TYPE );

  }

  return( WLZ_ERR_NONE );
}
