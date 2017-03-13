#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBackground_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzBackground.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
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
* \brief	Functions to get and set the background value of objects.
* \ingroup	WlzValuesUtils
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	Sets the backgound value of an image object.
* \param	obj			Given object in which to set the
*					background value.
* \param	bgd			Required background value.
*/
WlzErrorNum WlzSetBackground(WlzObject	*obj,
			     WlzPixelV	bgd)
{
  WlzPlaneDomain	*planedmn;
  WlzVoxelValues	*voxtab;
  WlzObject		*obj1;
  int			i;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  if( obj == NULL ){
    return( WLZ_ERR_OBJECT_NULL );
  }

  switch( obj->type ){

  case WLZ_2D_DOMAINOBJ:
  case WLZ_3D_DOMAINOBJ:
    break;

  case WLZ_TRANS_OBJ:
    return( WlzSetBackground(obj->values.obj, bgd) );

  case WLZ_EMPTY_OBJ:
    return( WLZ_ERR_NONE );

  default:
    return( WLZ_ERR_OBJECT_TYPE );

  }    

  if( obj->domain.core == NULL ){
    return( WLZ_ERR_DOMAIN_NULL );
  }
  if( obj->values.core == NULL ){
    return( WLZ_ERR_NONE );
  }

  /* check the background type */
  switch( bgd.type ){

  case WLZ_GREY_INT:
  case WLZ_GREY_SHORT:
  case WLZ_GREY_UBYTE:
  case WLZ_GREY_FLOAT:
  case WLZ_GREY_DOUBLE:
  case WLZ_GREY_RGBA:
    break;

  default:
    return( WLZ_ERR_GREY_TYPE );

  }
   
  switch( obj->type ){

  case WLZ_2D_DOMAINOBJ:
    (void )WlzValueConvertPixel(
				&bgd, bgd,
				WlzGreyTableTypeToGreyType(obj->values.core->type, &errNum));
    switch( WlzGreyTableTypeToTableType(obj->values.core->type, &errNum) ){

    case WLZ_GREY_TAB_RAGR:
      obj->values.v->bckgrnd = bgd;
      break;

    case WLZ_GREY_TAB_RECT:
      obj->values.r->bckgrnd = bgd;
      break;

    case WLZ_GREY_TAB_INTL:
      obj->values.i->bckgrnd = bgd;
      break;

    default:
      return( errNum );
    }
    break;

  case WLZ_3D_DOMAINOBJ:
    planedmn = obj->domain.p;
    if( planedmn->type != WLZ_PLANEDOMAIN_DOMAIN ){
      return( WLZ_ERR_PLANEDOMAIN_TYPE );
    }

    voxtab = obj->values.vox;
    if( voxtab->type != WLZ_VOXELVALUETABLE_GREY ){
      return( WLZ_ERR_VOXELVALUES_TYPE );
    }

    for(i=0; i <= planedmn->lastpl - planedmn->plane1; i++){
      if( planedmn->domains[i].core == NULL
	  || voxtab->values[i].core == NULL){
	continue;
      }
      obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, planedmn->domains[i],
			 voxtab->values[i], NULL, NULL, &errNum);
      WlzSetBackground( obj1, bgd );

      WlzFreeObj( obj1 );
    }
    voxtab->bckgrnd = bgd;
    break;

  default:
    errNum = WLZ_ERR_OBJECT_TYPE;
    break;
  }

  return( errNum );
}

/*!
* \return	The background value. If the returned pixel value type is
* 		WLZ_GREY_ERROR then an error has occurred.
* \ingroup	WlzValuesUtils
* \brief	Gets the background value of the given object.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPixelV WlzGetBackground(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzPlaneDomain	*planedmn;
  WlzPixelV		bgd;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* set up an invalid background return */
  bgd.type = WLZ_GREY_ERROR;
  bgd.v.inv = 0;
  
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }


  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.i == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }

      if( obj->values.v == NULL ){
	bgd.type = WLZ_GREY_INT;
	bgd.v.inv = 0;
	break;
      }
      
      switch( WlzGreyTableTypeToTableType(obj->values.core->type,
					  NULL) ){

      case WLZ_GREY_TAB_RAGR:
	bgd = obj->values.v->bckgrnd;
	break;

      case WLZ_GREY_TAB_RECT:
	bgd = obj->values.r->bckgrnd;
	break;

      case WLZ_GREY_TAB_INTL:
	bgd = obj->values.i->bckgrnd;
	break;

      case WLZ_GREY_TAB_TILED:
        bgd = obj->values.t->bckgrnd;
	break;

      default:
	errNum = WLZ_ERR_VALUES_TYPE;
	break;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.p == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }

      if( obj->values.core == NULL ){
	bgd.type = WLZ_GREY_INT;
	bgd.v.inv = 0;
	errNum = WLZ_ERR_NONE;
	break;
      }

      planedmn = obj->domain.p;
      if( planedmn->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
	break;
      }

      if(WlzGreyTableIsTiled(obj->values.core->type))
      {
        bgd = obj->values.t->bckgrnd;
      }
      else if(obj->values.core->type == WLZ_VOXELVALUETABLE_GREY)
      {
        bgd = obj->values.vox->bckgrnd;
      }
      else
      {
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      break;

    case WLZ_TRANS_OBJ:
      return( WlzGetBackground(obj->values.obj, dstErr) );

    case WLZ_EMPTY_OBJ:
      bgd.type = WLZ_GREY_INT;
      bgd.v.inv = 0;
      errNum = WLZ_ERR_NONE;
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return bgd;
}

/*!
* \return	New Woolz object or NULL on error. 
* \ingroup	WlzValuesUtils
* \brief	Creates a new Woolz object in which the background value
* 		is as given. The given object must be either a 2D or 3D
* 		domain object with a valid domain, however it's values can
* 		be NULL. If the values are NULL then a new value table is
* 		created with the same grey type as the given background
* 		value, or if non-NULL then a new value table is created
* 		while sharing the values of the original object's value
* 		table when this is possible. The given object remains
* 		unchanged apart from linkcounts.
* \param	gObj			Given object.
* \param	bgdV			Given background value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzSetBackGroundNewObj(
				  WlzObject *gObj,
				  WlzPixelV bgdV,
				  WlzErrorNum *dstErr)
{
  WlzObjectType tType;
  WlzObject	*rObj = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->values.core == NULL)
  {
    /* Given object has NULL values so create a new value table. */
    tType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, bgdV.type, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      rObj = WlzNewObjectValues(gObj, tType, bgdV, 1, bgdV, &errNum);
    }
  }
  else
  {
    /* Given object has values so create new value table using the same
     * values. */
    rObj = WlzNewObjectValueTable(gObj, bgdV, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}
