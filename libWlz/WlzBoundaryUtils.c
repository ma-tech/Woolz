#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBoundaryUtils_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzBoundaryUtils.c
* \author       Bill Hill, Richard Baldock
* \date         September 2002
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
* \brief	Simple functions that operate on boundary lists.
* \ingroup	WlzBoundary
*/

#include <Wlz.h>

static int			WlzBoundPolyCountFn(
				  WlzBoundList *bnd);
static void			WlzBoundObjToPolyFillArray(
				  WlzBoundList *bnd,
				  WlzPolygonDomain **polyAry,
				  int *idx);

/*!
* \return	Woolz error code.
* \ingroup	WlzBoundary
* \brief	decomposes a boundary into it's component polygons.
* \param	bndObj		Given boundary.
* \param	dstNumObjs	Destination pointer for the number of polygons.
* \param	dstObjArray	Destination pointer for the array of polygons.
*/
WlzErrorNum WlzBoundaryToPolyObjArray(
  WlzObject	*bndObj,
  int		*dstNumObjs,
  WlzObject	***dstObjArray)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzDomain	domain;
  WlzValues	values;
  WlzObject	*obj, **objs;
  WlzPolygonDomain	**polyArray;
  int 		i, numPolys;

  /* check inputs */
  if( bndObj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((dstNumObjs == NULL) || (dstObjArray == NULL)){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else {
    /* generate array of poly domains */
    errNum = WlzBoundObjToPolyDomArray(bndObj, &numPolys, &polyArray);
  }

  /* convert to polygon objects */
  if( errNum == WLZ_ERR_NONE ){
    if((objs = (WlzObject **) AlcMalloc(sizeof(WlzObject *)*numPolys)) == NULL){
      errNum = WLZ_ERR_MEM_ALLOC;
      for(i=0; i < numPolys; i++){
	WlzFreePolyDmn(polyArray[i]);
      }
      AlcFree(polyArray);
      numPolys = 0;
    }
    else {
      for(i=0; i < numPolys; i++){
	domain.poly = polyArray[i];
	values.core = NULL;
	obj = WlzMakeMain(WLZ_2D_POLYGON, domain, values,
			  NULL, NULL, &errNum);
	objs[i] = WlzAssignObject(obj, NULL);
	WlzFreePolyDmn(polyArray[i]);
      }
      AlcFree(polyArray);
    }
  }
  
  *dstNumObjs = numPolys;
  *dstObjArray = objs;
  return errNum;
}

/*!
* \return	Array of polygon domains.
* \ingroup	WlzBoundary
* \brief	Given a boundary list object returns a simple array of
*		polygon domains.
* \param	bndObj			Given boundary list object.
* \param	dstArySz		Destination ptr for array size.
* \param	dstPolyAry		Destination ptr for the array.
*/
WlzErrorNum	WlzBoundObjToPolyDomArray(WlzObject *bndObj, int *dstArySz,
					  WlzPolygonDomain ***dstPolyAry)
{
  int		idx,
  		polyCnt;
  WlzPolygonDomain **polyAry = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bndObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(bndObj->type != WLZ_BOUNDLIST)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(bndObj->domain.b == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstArySz == NULL) || (dstPolyAry == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  /* Count number of polydomains in the boundary list. */
  polyCnt = WlzBoundPolyCount(bndObj->domain.b, &errNum);
  /* Allocate array for polygon domain pointers. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((polyAry = (WlzPolygonDomain **)AlcMalloc(sizeof(WlzPolygonDomain *) *
    						 polyCnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Fill in the Array. */
  if(errNum == WLZ_ERR_NONE)
  {
    idx = 0;
    WlzBoundObjToPolyFillArray(bndObj->domain.b, polyAry, &idx);
    *dstArySz = polyCnt;
    *dstPolyAry = polyAry;
  }
  return(errNum);
}

/*!
* \return	The number of polygon domains in the given boundary list.
* \ingroup	WlzBoundary
* \brief	Count the number of polygon domains in a boundary list.
* \param	bnd			Given boundary list.
* \param	dstErr			Destination ptr for error, may be NULL.
*/
int		WlzBoundPolyCount(WlzBoundList *bnd, WlzErrorNum *dstErr)
{
  int		cnt = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bnd == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    cnt = WlzBoundPolyCountFn(bnd);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cnt);
}

/*!
* \return	The number of polygon domains in the given boundary list
*		from the given entry point downwards..
* \ingroup	WlzBoundary
* \brief	Recursive function which actualy does the counting of the
* 		number of polygon domains in the given boundary list and
*		its children.
* \param	bnd			Given boundary list.
*/
static int	WlzBoundPolyCountFn(WlzBoundList *bnd)
{
  int		cnt = 0;

  if(bnd)
  {
    do
    {
      ++cnt;
      cnt += WlzBoundPolyCountFn(bnd->down);
      bnd = bnd->next;
    } while(bnd != NULL);
  }
  return(cnt);
}

/*!
* \return
* \ingroup	WlzBoundary
* \brief	Recursive function which counts the total number of
*		vertices in a boundary list.
* \param	bnd			Given boundary list.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzBoundVtxCount(WlzBoundList *bnd, WlzErrorNum *dstErr)
{
  int		cnt = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bnd == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    if(bnd->poly)
    {
      cnt += bnd->poly->nvertices;
    }
    if(bnd->next)
    {
      cnt += WlzBoundVtxCount(bnd->next, &errNum);
    }
    if(bnd->down)
    {
      cnt += WlzBoundVtxCount(bnd->down, &errNum);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cnt);
}

/*!
* \return
* \ingroup	WlzBoundary
* \brief	Sets polygon domains from a boundary list to the given
*		array, which must be large enough. This is a recursive
*		function.
* \param	bnd			Given boundary list.
* \param	polyAry			The polygon domains array.
* \param	idx			Index incremented for each assignment.
*/
static void	WlzBoundObjToPolyFillArray(WlzBoundList *bnd,
					   WlzPolygonDomain **polyAry,
					   int *idx)
{
  if(bnd)
  {
    do
    {
      *(polyAry + *idx) = WlzAssignPolygonDomain(bnd->poly, NULL);
      ++*idx;
      WlzBoundObjToPolyFillArray(bnd->down, polyAry, idx);
      bnd = bnd->next;
    } while(bnd != NULL);
  }
}
