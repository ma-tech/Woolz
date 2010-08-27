#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFTxt_c[] = "MRC HGU $Id:";
#endif
#endif
/*!
* \file         WlzExtFFTxt.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Apr 28 14:53:55 2010
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par Copyright:
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
* \brief	Functions for writing Woolz objects to 
*		a simple text '.txt' data format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>


WlzErrorNum WlzEffWriteObjTxt(
  FILE 		*fP,
  WlzObject 	*obj,
  char		*params)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		width, height;
  int		i, j;
  WlzObject	*rectObj=NULL;
  WlzGreyType 	gType;
  WlzGreyP	gValuesP;

  /* check input */
  if( fP == NULL ){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core ){
	if( obj->values.core == NULL ){
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else {
	  WlzIBox2	cutBox;
	  
	  cutBox.xMin = obj->domain.i->kol1;
	  cutBox.yMin = obj->domain.i->line1;
	  cutBox.xMax = obj->domain.i->lastkl;
	  cutBox.yMax = obj->domain.i->lastln;
	  gType = WlzGreyTypeFromObj(obj, &errNum);
	  if((rectObj = WlzCutObjToBox2D(obj, cutBox, gType,
					 0, 0.0, 0.0, &errNum)) != NULL){
	    width = rectObj->domain.i->lastkl - rectObj->domain.i->kol1 + 1;
	    height = rectObj->domain.i->lastln - rectObj->domain.i->line1 + 1;
	    gValuesP = rectObj->values.r->values;
	  }
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_TRANS_OBJ:
      errNum = WlzEffWriteObjTxt(fP, obj->values.obj, params);
      break;

    case WLZ_EMPTY_OBJ:
      return WLZ_ERR_NONE;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* now write each row with comma-separated values */
  if( errNum == WLZ_ERR_NONE ){
    for(j=0; j < height; j++ ){
      for(i=0; i < (width-1); i++){
	switch( gType ){
	case WLZ_GREY_UBYTE:
	  fprintf(fP, "%d, ", *gValuesP.ubp);
	  gValuesP.ubp++;
	  break;
	case WLZ_GREY_SHORT:
	  fprintf(fP, "%d, ", *gValuesP.shp);
	  gValuesP.shp++;
	  break;
	case WLZ_GREY_INT:
	  fprintf(fP, "%d, ", *gValuesP.inp);
	  gValuesP.inp++;
	  break;
	default:
	  break;
	}
      }
      switch( gType ){
      case WLZ_GREY_UBYTE:
	fprintf(fP, "%d\n", *gValuesP.ubp);
	gValuesP.ubp++;
	break;
      case WLZ_GREY_SHORT:
	fprintf(fP, "%d\n", *gValuesP.shp);
	gValuesP.shp++;
	break;
      case WLZ_GREY_INT:
	fprintf(fP, "%d\n", *gValuesP.inp);
	gValuesP.inp++;
	break;
      default:
	break;
      }
    }
    WlzFreeObj(rectObj);
  }

  return(errNum);
}
