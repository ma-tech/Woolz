#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConstruct3D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzConstruct3D.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Functions to construct 3D domain objects from 2D domain
* 		objects.
* \ingroup	WlzAllocation
*/

#include <stdio.h>
#include <Wlz.h>

/*!
* \return	New 3D object.
* \ingroup	WlzAllocation
* \brief	Constructs a 3D domain object from 2D domain objects read
*		from the given files. Each file is read in turn and added
*		to the 3D object. An empty plane can be specified by
*		setting the file string to NULL. Either all or none of
*		the 2D objects must have values. When the 2D objects
*		have values then the background value of the first 2D
*		object is set to be the background value of the 3D object.
* \param	nFileStr		Number of file strings.
* \param	fileStr			File strings.
* \param	plane1			The plane coordinate of the first
*					2D object.
* \param	xSz			Column voxel size.
* \param	ySz			Line voxel size.
* \param	zSz			Plane voxel size.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzConstruct3DObjFromFile(int nFileStr, char **fileStr,
					  int plane1,
					  float xSz, float ySz, float zSz,
					  WlzErrorNum *dstErr)
{
  int		idx,
  		lastpl;
  WlzDomain	dom3D;
  WlzValues	val3D;
  WlzObject	*obj2D = NULL,
  		*obj3D = NULL;
  WlzPixelV	bgd;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom3D.core = NULL;
  val3D.core = NULL;
  lastpl = plane1 + nFileStr - 1;
  if((nFileStr <= 0) || (fileStr == NULL) || (*fileStr == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if((fP = fopen(*fileStr, "r")) == NULL)
    {
      errNum = WLZ_ERR_READ_EOF;
    }
    else
    {
      obj2D = WlzReadObj(fP, &errNum);
      (void )fclose(fP);
      fP = NULL;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj2D->type)
    {
      case WLZ_EMPTY_OBJ:
	break;
      case WLZ_2D_DOMAINOBJ:
	if(obj2D->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Make a plane domain, set column and line bounds later. */
  if(errNum == WLZ_ERR_NONE)
  {
    dom3D.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
    				 plane1, lastpl,
				 0, 1, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom3D.p->voxel_size[0] = xSz;
    dom3D.p->voxel_size[1] = ySz;
    dom3D.p->voxel_size[2] = zSz;
  }
  /* Make a voxel value table. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(obj2D->values.core)
    {
      bgd = WlzGetBackground(obj2D, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	val3D.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					plane1, lastpl, bgd, NULL, &errNum);
      }
    }				    
  }
  idx = 0;
  while((errNum == WLZ_ERR_NONE) && (idx < nFileStr))
  {
    if(obj2D)
    {
      switch(obj2D->type)
      {
	case WLZ_EMPTY_OBJ:
	  break;
	case WLZ_2D_DOMAINOBJ:
	  if(obj2D->domain.core == NULL)
	  {
	    errNum = WLZ_ERR_DOMAIN_NULL;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	*(dom3D.p->domains + idx) = WlzAssignDomain(obj2D->domain, NULL);
	if(val3D.core)
	{
	  if((obj2D->domain.core != NULL) && (obj2D->values.core == NULL))
	  {
	    errNum = WLZ_ERR_VALUES_NULL;
	  }
	  else
	  {
	    *(val3D.vox->values + idx) = WlzAssignValues(obj2D->values, NULL);
	  }
	}
      }
      WlzFreeObj(obj2D);
      obj2D = NULL;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      ++idx;
      if((idx < nFileStr) && *(fileStr + idx))
      {
	if((fP = fopen(*(fileStr + idx), "r")) == NULL)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
	else
	{
	  obj2D = WlzReadObj(fP, &errNum);
	  (void )fclose(fP);
	  fP = NULL;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(dom3D.p, val3D.vox);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj3D = WlzMakeMain(WLZ_3D_DOMAINOBJ, dom3D, val3D, NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dom3D.core)
    {
      (void )WlzFreeDomain(dom3D);
    }
    if(val3D.core)
    {
      (void )WlzFreeValues(val3D);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj3D);
}

/*!
* \return	New 3D object.
* \ingroup	WlzAllocation
* \brief	Constructs a 3D domain object from 2D domain objects. Each
*		2D object is assigned in turn to the 3D object no domains
*		or values are copied. An empty can be specified by
*		setting the 2D object to NULL. Either all or none of
*		the 2D objects must have values. When the 2D objects
*		have values then the background value of the first 2D
*		object is set to be the background value of the 3D object.
* \param	nObjs			Number of objects.
* \param	objs			The 2D objects, the first of which
*					MUST not be NULL.
* \param	plane1			The plane coordinate of the first
*					2D object.
* \param	xSz			Column voxel size.
* \param	ySz			Line voxel size.
* \param	zSz			Plane voxel size.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzConstruct3DObjFromObj(int nObjs, WlzObject **objs,
					 int plane1,
					 float xSz, float ySz, float zSz,
					 WlzErrorNum *dstErr)
{
  int		idx,
  		lastpl;
  WlzDomain	dom3D;
  WlzValues	val3D;
  WlzObject	*obj2D,
  		*obj3D = NULL;
  WlzPixelV	bgd;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom3D.core = NULL;
  val3D.core = NULL;
  lastpl = plane1 + nObjs - 1;
  if((nObjs <= 0) || (objs == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((obj2D = *objs) == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(obj2D->type)
    {
      case WLZ_EMPTY_OBJ:
	break;
      case WLZ_2D_DOMAINOBJ:
	if(obj2D->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Make a plane domain, set column and line bounds later. */
  if(errNum == WLZ_ERR_NONE)
  {
    dom3D.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
    				 plane1, lastpl,
				 0, 1, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom3D.p->voxel_size[0] = xSz;
    dom3D.p->voxel_size[1] = ySz;
    dom3D.p->voxel_size[2] = zSz;
  }
  /* Make a voxel value table. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(obj2D->values.core)
    {
      bgd = WlzGetBackground(obj2D, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	val3D.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					plane1, lastpl, bgd, NULL, &errNum);
      }
    }				    
  }
  for(idx = 0; (errNum == WLZ_ERR_NONE) && (idx < nObjs); ++idx)
  {
    obj2D = *(objs + idx);
    if(obj2D)
    {
      switch(obj2D->type)
      {
	case WLZ_EMPTY_OBJ:
	  break;
	case WLZ_2D_DOMAINOBJ:
	  if(obj2D->domain.core == NULL)
	  {
	    errNum = WLZ_ERR_DOMAIN_NULL;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	*(dom3D.p->domains + idx) = WlzAssignDomain(obj2D->domain, NULL);
	if(val3D.core)
	{
	  if((obj2D->domain.core != NULL) && (obj2D->values.core == NULL))
	  {
	    errNum = WLZ_ERR_VALUES_NULL;
	  }
	  else
	  {
	    *(val3D.vox->values + idx) = WlzAssignValues(obj2D->values, NULL);
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(dom3D.p, val3D.vox);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj3D = WlzMakeMain(WLZ_3D_DOMAINOBJ, dom3D, val3D, NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dom3D.core)
    {
      (void )WlzFreeDomain(dom3D);
    }
    if(val3D.core)
    {
      (void )WlzFreeValues(val3D);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj3D);
}
