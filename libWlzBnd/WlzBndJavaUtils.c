#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBndJavaUtils_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzBnd/WlzBndJavaUtils.c
* \author       Nick Burton, Bill Hill
* \date         August 2003
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
* \brief	Binding functions to enable Java Woolz to access the
*		Woolz data structures.
* \ingroup	LibWlzBnd
*/
#include <Wlz.h>

/*!
* \return	Woolz object type as an int.
* \ingroup	LibWlzBnd
* \brief	Gets the Woolz object type.
* \param	obj		Given Woolz object.
*/
int WlzBndObjGetType(WlzObject *obj)
{
  return(obj->type);
}

/*!
* \return	First line of 2D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 2D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the first line of a Woolz 2D domain object.
* \param	obj
*/
int WlzBndGetIntervalDomainLine1(WlzObject *obj)
{
   int	line1 = -1;

   if(obj)
   {
     line1 = (obj->type == WLZ_2D_DOMAINOBJ)? obj->domain.i->line1: -2;
   }
   return(line1);
}

/*!
* \return	Last line of 2D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 2D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the last line of a Woolz 2D domain object.
* \param	obj
*/
int WlzBndGetIntervalDomainLastLn(WlzObject *obj)
{
   int	lastln = -1;

   if(obj)
   {
     lastln = (obj->type == WLZ_2D_DOMAINOBJ)? obj->domain.i->lastln: -2;
   }
   return(lastln);
}

/*!
* \return	First column of 2D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 2D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the first column of a Woolz 2D domain object.
* \param	obj
*/
int WlzBndGetIntervalDomainKol1(WlzObject *obj)
{
   int	kol1 = -1;

   if(obj)
   {
     kol1 = (obj->type == WLZ_2D_DOMAINOBJ)? obj->domain.i->kol1: -2;
   }
   return(kol1);
}

/*!
* \return	Last column of 2D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 2D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the last column of a Woolz 2D domain object.
* \param	obj
*/
int WlzBndGetIntervalDomainLastKl(WlzObject *obj)
{
   int	lastkl = -1;

   if(obj)
   {
     lastkl = (obj->type == WLZ_2D_DOMAINOBJ)? obj->domain.i->lastkl: -2;
   }
   return(lastkl);
}

/*!
* \return	First line of 3D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 3D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the first line of a Woolz 3D domain object.
* \param	obj
*/
int WlzBndGetPlaneDomainLine1(WlzObject *obj)
{
   int	line1 = -1;

   if(obj)
   {
     line1 = (obj->type == WLZ_3D_DOMAINOBJ)? obj->domain.p->line1: -2;
   }
   return(line1);
}

/*!
* \return	Last line of 3D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 3D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the last line of a Woolz 3D domain object.
* \param	obj
*/
int WlzBndGetPlaneDomainLastLn(WlzObject *obj)
{
   int	lastln = -1;

   if(obj)
   {
     lastln = (obj->type == WLZ_3D_DOMAINOBJ)? obj->domain.p->lastln: -2;
   }
   return(lastln);
}

/*!
* \return	First column of 3D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 3D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the first column of a Woolz 3D domain object.
* \param	obj
*/
int WlzBndGetPlaneDomainKol1(WlzObject *obj)
{
   int	kol1 = -1;

   if(obj)
   {
     kol1 = (obj->type == WLZ_3D_DOMAINOBJ)? obj->domain.p->kol1: -2;
   }
   return(kol1);
}

/*!
* \return	Last column of 3D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 3D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the last column of a Woolz 3D domain object.
* \param	obj
*/
int WlzBndGetPlaneDomainLastKl(WlzObject *obj)
{
   int	lastkl = -1;

   if(obj)
   {
     lastkl = (obj->type == WLZ_3D_DOMAINOBJ)? obj->domain.p->lastkl: -2;
   }
   return(lastkl);
}

/*!
* \return	First plane of 3D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 3D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the first plane of a Woolz 3D domain object.
* \param	obj
*/
int WlzBndGetPlaneDomainPlane1(WlzObject *obj)
{
   int	plane1 = -1;

   if(obj)
   {
     plane1 = (obj->type == WLZ_3D_DOMAINOBJ)? obj->domain.p->plane1: -2;
   }
   return(plane1);
}

/*!
* \return	Last plane of 3D domain object if object valid,
* 		-1 if object NULL or -2 if object not a 3D
* 		domain object.
* \ingroup	LibWlzBnd
* \brief	Gets the last plane of a Woolz 3D domain object.
* \param	obj
*/
int WlzBndGetPlaneDomainLastPl(WlzObject *obj)
{
   int	lastpl = -1;

   if(obj)
   {
     lastpl = (obj->type == WLZ_3D_DOMAINOBJ)? obj->domain.p->lastpl: -2;
   }
   return(lastpl);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	LibWlzBnd
* \brief	Creates a new Woolz object with the domain of the first
* 		object and the values of the second.
* \param	domainObj		Object with domain to use.
* \param	valuesObj		Object with values to use.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzBndSetGreyValues(WlzObject *domainObj,
			WlzObject *valuesObj,
			WlzErrorNum *dstErr) 
{

  WlzObject *retObj = NULL;
  WlzErrorNum  errNum = WLZ_ERR_NONE;

  if(domainObj && (domainObj->type == WLZ_2D_DOMAINOBJ) &&
     valuesObj && (valuesObj->type == WLZ_2D_DOMAINOBJ))
  {
    retObj = WlzMakeMain(domainObj->type, domainObj->domain, valuesObj->values,
                         NULL, NULL, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(retObj);
}

/*!
* \return 	Affine transform of the given view struct.
* \ingroup	LibWlzBnd
* \brief	Gets the affine transform of the given view struct.
* \param	viewStr			Given view struct.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzAffineTransform *WlzBndGetAffineTransform(WlzThreeDViewStruct *viewStr,
    				WlzErrorNum *dstErr)
{

  WlzAffineTransform *retObj = NULL;
  WlzErrorNum  errNum = WLZ_ERR_NONE;

  if(viewStr)
  {
    retObj = viewStr->trans;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(retObj);
}

/*!
* \ingroup	LibWlzBnd
* \brief	Gets the affine transform matrix with size.
* \param	dstSizeArrayDat		Destination pointer for array data.
* \param	dstArrayDat		Destination pointer for array size,
* 					always 4 x 4.
* \param	trans			Given affine transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
void 		WlzBndGetTransformMatrix(WlzIVertex2 *dstSizeArrayDat,
				double  ***dstArrayDat,
				WlzAffineTransform *trans,
				WlzErrorNum *dstErr)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzIVertex2 size;

  size.vtX = 4;
  size.vtY = 4;
  if(dstArrayDat == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    *dstArrayDat = trans->mat;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
}
